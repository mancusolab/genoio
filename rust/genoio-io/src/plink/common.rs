// pattern: Functional Core
//! Shared PLINK metadata normalization.
//!
//! PLINK1 and PLINK2 use different missing-value tokens in companion text
//! files. The parsers call these helpers before constructing core metadata
//! records.

pub(crate) const PLINK1_MISSING_VALUES: &[&str] = &["0"];
pub(crate) const PLINK2_MISSING_VALUES: &[&str] = &["0", ".", "NA"];

pub(crate) fn optional_plink_value(value: &str, missing_values: &[&str]) -> Option<String> {
    if missing_values.contains(&value) {
        None
    } else {
        Some(value.to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn optional_plink_value_uses_format_specific_missing_tokens() {
        for missing_value in PLINK1_MISSING_VALUES {
            assert_eq!(
                optional_plink_value(missing_value, PLINK1_MISSING_VALUES),
                None
            );
        }
        assert_eq!(
            optional_plink_value(".", PLINK1_MISSING_VALUES),
            Some(".".to_string())
        );
        assert_eq!(
            optional_plink_value("NA", PLINK1_MISSING_VALUES),
            Some("NA".to_string())
        );

        for missing_value in PLINK2_MISSING_VALUES {
            assert_eq!(
                optional_plink_value(missing_value, PLINK2_MISSING_VALUES),
                None
            );
        }
        assert_eq!(
            optional_plink_value("sample-1", PLINK2_MISSING_VALUES),
            Some("sample-1".to_string())
        );
    }
}
