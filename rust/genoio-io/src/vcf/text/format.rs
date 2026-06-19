//! Shared byte scanners for VCF FORMAT/sample columns.
//!
//! The text backend decoders use these helpers to find a FORMAT key once and then
//! visit only selected sample columns without materializing intermediate fields.

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum FormatScanError<E> {
    /// The FORMAT/sample structure was malformed before a decoder saw a token.
    Scan(&'static str),
    /// The caller rejected a token after FORMAT/sample navigation succeeded.
    Emit(E),
}

pub(super) fn format_key_index(sample_fields: &[u8], key: &[u8]) -> Option<usize> {
    let format_end = sample_fields.iter().position(|&b| b == b'\t')?;
    sample_fields[..format_end]
        .split(|&b| b == b':')
        .position(|candidate| candidate == key)
}

pub(super) fn scan_selected_format_tokens<E>(
    sample_fields: &[u8],
    key_index: usize,
    source_indices: &[usize],
    missing_value_error: &'static str,
    emit: &mut impl FnMut(&[u8]) -> std::result::Result<(), E>,
) -> std::result::Result<(), FormatScanError<E>> {
    let Some(format_end) = sample_fields.iter().position(|&b| b == b'\t') else {
        return Err(FormatScanError::Scan(
            "record has FORMAT but no sample columns",
        ));
    };
    let mut selected_index = 0_usize;
    let mut sample_index = 0_usize;
    let mut field_start = format_end + 1;

    while selected_index < source_indices.len() {
        let target_index = source_indices[selected_index];
        if field_start > sample_fields.len() {
            return Err(FormatScanError::Scan(
                "selected sample index is outside the record",
            ));
        }
        let field_end = next_delimiter(sample_fields, field_start, b'\t');
        if sample_index == target_index {
            let token = nth_colon_field(&sample_fields[field_start..field_end], key_index)
                .ok_or(FormatScanError::Scan(missing_value_error))?;
            emit(token).map_err(FormatScanError::Emit)?;
            selected_index += 1;
        }
        sample_index += 1;
        if field_end == sample_fields.len() {
            field_start = sample_fields.len() + 1;
        } else {
            field_start = field_end + 1;
        }
    }

    Ok(())
}

pub(super) fn nth_colon_field(sample: &[u8], index: usize) -> Option<&[u8]> {
    let mut field_start = 0_usize;
    for field_index in 0..=index {
        let field_end = next_delimiter(sample, field_start, b':');
        if field_index == index {
            return Some(&sample[field_start..field_end]);
        }
        if field_end == sample.len() {
            return None;
        }
        field_start = field_end + 1;
    }
    None
}

pub(super) fn next_delimiter(buf: &[u8], start: usize, delimiter: u8) -> usize {
    buf[start..]
        .iter()
        .position(|&b| b == delimiter)
        .map_or(buf.len(), |offset| start + offset)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn format_key_index_finds_requested_key() {
        assert_eq!(format_key_index(b"GT:DS\t0/1:0.4", b"GT"), Some(0));
        assert_eq!(format_key_index(b"DP:GT:GQ\t9:1/1:99", b"GT"), Some(1));
        assert_eq!(format_key_index(b"DP:GQ\t9:99", b"GT"), None);
    }

    #[test]
    fn scan_selected_format_tokens_visits_source_order_subset() {
        let mut tokens = Vec::new();
        scan_selected_format_tokens(
            b"DP:DS\t9:0.0\t8:1.0\t7:2.0\t6:.",
            1,
            &[1, 3],
            "sample is missing DS value",
            &mut |token| {
                tokens.push(token.to_vec());
                Ok::<(), &'static str>(())
            },
        )
        .expect("selected tokens should scan");

        assert_eq!(tokens, vec![b"1.0".to_vec(), b".".to_vec()]);
    }

    #[test]
    fn scan_selected_format_tokens_reports_missing_key_value() {
        let error = scan_selected_format_tokens(
            b"GT:DS\t0/1",
            1,
            &[0],
            "sample is missing DS value",
            &mut |_token| Ok::<(), &'static str>(()),
        )
        .expect_err("short sample field should fail");

        assert_eq!(error, FormatScanError::Scan("sample is missing DS value"));
    }
}
