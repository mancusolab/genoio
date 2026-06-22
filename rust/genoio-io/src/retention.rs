// pattern: Functional Core
//! Retained-variant window state for metadata and genotype filters.
//!
//! The state counts variants that survive filtering, not raw source rows. This
//! keeps source-window reads and retained-window reads separate in the callers.

use genoio_core::{DenseDiagnostics, PartialFilterDecision, VariantWindow};

/// Action after evaluating metadata-only predicates for one source variant.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum MetadataRetentionAction {
    Include,
    DecodeGenotypes,
    Skip,
    Stop,
}

/// Action after a variant is known to survive or fail genotype-stat filters.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum RetentionAction {
    Include,
    Skip,
    Stop,
}

/// Tracks retained-variant index for optional retained-output windows.
///
/// Metadata rejects do not advance the retained index. Metadata accepts and
/// genotype accepts do, because those variants are retained before the optional
/// window is applied.
#[derive(Debug, Clone, Copy)]
pub(crate) struct RetainedVariantState {
    window: Option<VariantWindow>,
    retained_index: usize,
}

impl RetainedVariantState {
    pub(crate) fn new(window: Option<VariantWindow>) -> Self {
        Self {
            window,
            retained_index: 0,
        }
    }

    pub(crate) fn metadata_decision(
        &mut self,
        partial_decision: PartialFilterDecision,
        diagnostics: &mut DenseDiagnostics,
    ) -> MetadataRetentionAction {
        diagnostics.candidate_variants += 1;
        match partial_decision {
            PartialFilterDecision::Reject => {
                diagnostics.dropped_metadata_variants += 1;
                MetadataRetentionAction::Skip
            }
            PartialFilterDecision::Accept => match self.consume_retained_variant() {
                RetentionAction::Include => MetadataRetentionAction::Include,
                RetentionAction::Skip => MetadataRetentionAction::Skip,
                RetentionAction::Stop => MetadataRetentionAction::Stop,
            },
            PartialFilterDecision::NeedGenotypes => MetadataRetentionAction::DecodeGenotypes,
        }
    }

    pub(crate) fn genotype_decision(
        &mut self,
        retain_variant: bool,
        diagnostics: &mut DenseDiagnostics,
    ) -> RetentionAction {
        if !retain_variant {
            diagnostics.dropped_genotype_variants += 1;
            return RetentionAction::Skip;
        }
        self.consume_retained_variant()
    }

    pub(crate) fn window_is_satisfied(&self) -> bool {
        self.window
            .is_some_and(|window| window.is_past(self.retained_index))
    }

    fn consume_retained_variant(&mut self) -> RetentionAction {
        let include_in_window = self
            .window
            .is_none_or(|window| window.contains(self.retained_index));
        self.retained_index += 1;
        if include_in_window {
            RetentionAction::Include
        } else if self.window_is_satisfied() {
            RetentionAction::Stop
        } else {
            RetentionAction::Skip
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn metadata_reject_counts_drop_without_advancing_retained_index() {
        let mut state = RetainedVariantState::new(Some(VariantWindow { start: 0, len: 1 }));
        let mut diagnostics = DenseDiagnostics::default();

        assert_eq!(
            state.metadata_decision(PartialFilterDecision::Reject, &mut diagnostics),
            MetadataRetentionAction::Skip
        );
        assert_eq!(diagnostics.candidate_variants, 1);
        assert_eq!(diagnostics.dropped_metadata_variants, 1);
        assert!(!state.window_is_satisfied());
    }

    #[test]
    fn metadata_accept_skips_until_window_start_then_includes() {
        let mut state = RetainedVariantState::new(Some(VariantWindow { start: 1, len: 1 }));
        let mut diagnostics = DenseDiagnostics::default();

        assert_eq!(
            state.metadata_decision(PartialFilterDecision::Accept, &mut diagnostics),
            MetadataRetentionAction::Skip
        );
        assert_eq!(
            state.metadata_decision(PartialFilterDecision::Accept, &mut diagnostics),
            MetadataRetentionAction::Include
        );
        assert!(state.window_is_satisfied());
    }

    #[test]
    fn genotype_reject_does_not_advance_retained_index() {
        let mut state = RetainedVariantState::new(Some(VariantWindow { start: 0, len: 1 }));
        let mut diagnostics = DenseDiagnostics::default();

        assert_eq!(
            state.metadata_decision(PartialFilterDecision::NeedGenotypes, &mut diagnostics),
            MetadataRetentionAction::DecodeGenotypes
        );
        assert_eq!(
            state.genotype_decision(false, &mut diagnostics),
            RetentionAction::Skip
        );
        assert!(!state.window_is_satisfied());
        assert_eq!(
            state.genotype_decision(true, &mut diagnostics),
            RetentionAction::Include
        );
        assert!(state.window_is_satisfied());
        assert_eq!(diagnostics.dropped_genotype_variants, 1);
    }

    #[test]
    fn zero_length_window_stops_on_first_retained_variant() {
        let mut state = RetainedVariantState::new(Some(VariantWindow { start: 0, len: 0 }));
        let mut diagnostics = DenseDiagnostics::default();

        assert_eq!(
            state.metadata_decision(PartialFilterDecision::Accept, &mut diagnostics),
            MetadataRetentionAction::Stop
        );
    }
}
