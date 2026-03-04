import numpy as np
import pandas as pd
import pytest
from multiHGtest import from_time_to_event_to_survival_table, hg_test_dashboard
from multiHGtest import testHG_dashboard as _testHG_dashboard


class TestDashboard:
    """Test cases for hg_test_dashboard and the testHG_dashboard alias."""

    def test_dashboard_from_time_to_event(self):
        """Test full pipeline: time-to-event data -> survival table -> hg_test_dashboard."""
        data = pd.DataFrame({
            'time':  [1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6,
                      1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6],
            'event': [1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0,
                      1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1],
            'group': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                      2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
        })

        result = from_time_to_event_to_survival_table(data)
        Nt1 = result['at_risk_1']
        Nt2 = result['at_risk_2']
        Ot1 = result['observed_1']
        Ot2 = result['observed_2']

        df, stats = hg_test_dashboard(Nt1, Nt2, Ot1, Ot2)

        # DataFrame should have one row per time point
        assert len(df) == len(result['time'])

        # Check expected columns
        for col in ['at-risk1', 'at-risk2', 'events1', 'events2', 'pvalue', 'HCT']:
            assert col in df.columns, f"Missing column: {col}"

        # at-risk and event columns should match the inputs
        np.testing.assert_array_equal(df['at-risk1'].values, Nt1)
        np.testing.assert_array_equal(df['at-risk2'].values, Nt2)
        np.testing.assert_array_equal(df['events1'].values, Ot1)
        np.testing.assert_array_equal(df['events2'].values, Ot2)

        # P-values should be in [0, 1]
        assert np.all(df['pvalue'] >= 0) and np.all(df['pvalue'] <= 1)

        # HCT column should be boolean
        assert df['HCT'].dtype == bool

        # Stats dict should contain valid entries
        for key in ['hc', 'fisher', 'fisher_pval', 'minP']:
            assert key in stats, f"Missing stat: {key}"
            assert isinstance(stats[key], (int, float, np.floating))
            assert not np.isnan(stats[key])

        assert 0 <= stats['fisher_pval'] <= 1
        assert not np.isinf(stats['hc'])

    def test_dashboard_with_direct_arrays(self):
        """Test hg_test_dashboard with directly provided arrays."""
        Nt1 = np.array([100, 95, 90, 85])
        Nt2 = np.array([100, 92, 88, 82])
        Ot1 = np.array([5, 5, 5, 3])
        Ot2 = np.array([8, 4, 6, 3])

        df, stats = hg_test_dashboard(Nt1, Nt2, Ot1, Ot2)

        assert len(df) == 4
        assert 'pvalue' in df.columns
        assert 'HCT' in df.columns
        assert 'hc' in stats
        assert 'fisher' in stats
        assert 'fisher_pval' in stats
        assert 'minP' in stats

    def test_dashboard_alias(self):
        """Test that testHG_dashboard is an alias for hg_test_dashboard."""
        Nt1 = np.array([100, 95, 90, 85])
        Nt2 = np.array([100, 92, 88, 82])
        Ot1 = np.array([5, 5, 5, 3])
        Ot2 = np.array([8, 4, 6, 3])

        df1, stats1 = hg_test_dashboard(Nt1, Nt2, Ot1, Ot2)
        df2, stats2 = _testHG_dashboard(Nt1, Nt2, Ot1, Ot2)

        np.testing.assert_array_equal(df1['pvalue'].values, df2['pvalue'].values)
        assert stats1['hc'] == stats2['hc']
        assert stats1['fisher'] == stats2['fisher']

    def test_dashboard_pvals_alternative(self):
        """Test dashboard with different pvals_alternative options."""
        Nt1 = np.array([100, 95, 90, 85])
        Nt2 = np.array([100, 92, 88, 82])
        Ot1 = np.array([5, 5, 5, 3])
        Ot2 = np.array([8, 4, 6, 3])

        for alt in ['greater', 'less', 'two-sided']:
            df, stats = hg_test_dashboard(Nt1, Nt2, Ot1, Ot2, pvals_alternative=alt)
            assert np.all(df['pvalue'] >= 0) and np.all(df['pvalue'] <= 1)
            assert 0 <= stats['fisher_pval'] <= 1
