"""Reporting policy tests using known patient correlations and planning inputs."""
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import integrate, special, stats

import pase


def correlations(values, n=40):
    x = np.tile([-1., 1., -1., 1.], n // 4)
    e = np.tile([-1., -1., 1., 1.], n // 4)
    blocks = [np.column_stack([x, r*x + np.sqrt(1-r*r)*e]) for r in values]
    obs = pd.DataFrame({'donor': np.repeat([f'p{i}' for i in range(len(values))], n)},
                       index=[f'c{i}' for i in range(n*len(values))])
    return AnnData(np.vstack(blocks), obs=obs, var=pd.DataFrame(index=['A', 'B']))


def call(values, **kwargs):
    return pase.associate(correlations(values), 'A', 'B', 'donor', **kwargs)


def test_default_floor_withholds_five_and_retains_descriptions():
    res = call([-.6, -.2, .1, .2, .4], planning_sd=.2)
    s = res.summary().iloc[0]
    assert s.status == 'withheld_policy'
    assert s.reporting_status == 'descriptive_only'
    assert not s.inference_available
    assert s.min_patients == 6
    assert s.descriptive_mean_r == pytest.approx(-.02)
    assert s.descriptive_median_r == pytest.approx(.1)
    assert (s.n_positive, s.n_negative, s.n_zero) == (3, 2, 0)
    assert np.isfinite(res.patients.r).all()
    for key in ('estimate', 'ci_low', 'ci_high', 'p_value', 'standard_error',
                'df', 'tau2', 'planning_mde'):
        assert pd.isna(s[key]), key
    assert s.mde_status == 'withheld_policy'
    assert res.metadata['calibration_guaranteed'] is False


def test_six_patients_report_with_warning_and_keep_manual_t_result():
    values = [-.6, -.2, .1, .2, .4, .6]
    res = call(values)
    s = res.summary().iloc[0]
    se = np.std(values, ddof=1)/np.sqrt(6)
    assert s.status == 'reported' and s.inference_available
    assert s.reporting_status == 'inferential_with_warnings'
    assert s.estimate == pytest.approx(np.mean(values))
    assert s.ci_low == pytest.approx(np.mean(values)-stats.t.ppf(.975, 5)*se)
    assert 'small_cohort_calibration_warning' in s.warning_codes
    assert 'wide_interval' in s.warning_codes
    assert s.ci_width == pytest.approx(s.ci_high-s.ci_low)
    assert s.warning_codes == tuple(res.metadata['warning_codes'])


def test_method_specific_eligibility_controls_floor():
    values = [-.5, -.2, .1, .3, .6, 1.]
    equal = call(values)
    fisher = call(values, aggregate='re_hksj')
    assert equal.summary().iloc[0].inference_available
    assert fisher.summary().iloc[0].n_patients_eligible == 5
    assert fisher.summary().iloc[0].reporting_status == 'descriptive_only'
    assert fisher.patients.iloc[-1].reason == 'boundary_r'


def test_explicit_lower_floor_is_visible_override():
    res = call([-.3, .1, .4, .5], min_patients=2)
    s = res.summary().iloc[0]
    assert s.status == 'reported'
    assert 'below_default_min_patients' in s.warning_codes
    assert res.metadata['reporting_policy']['custom_min_patients']
    assert res.metadata['calibration_guaranteed'] is False


def test_warning_floor_is_not_a_calibration_certificate():
    s = call(np.linspace(-.3, .5, 20), max_ci_width=None)
    assert 'small_cohort_calibration_warning' not in s.summary().iloc[0].warning_codes
    assert s.metadata['calibration_guaranteed'] is False
    assert 'not a calibration guarantee' in s.metadata['validation_scope']


def test_width_threshold_warns_without_changing_inference():
    values = np.linspace(-.4, .6, 8)
    one = call(values, max_ci_width=None).summary().iloc[0]
    two = call(values, max_ci_width=.01).summary().iloc[0]
    assert 'wide_interval' not in one.warning_codes
    assert 'wide_interval' in two.warning_codes
    for k in ('estimate', 'ci_low', 'ci_high', 'p_value', 'standard_error'):
        assert one[k] == two[k]


def test_noncentral_t_mde_matches_independent_power_equation():
    s = call(np.linspace(-.3, .5, 20), planning_sd=.2,
             minimum_relevant_effect=.03).summary().iloc[0]
    critical = stats.t.ppf(.975, 19)
    ncp = s.planning_mde*np.sqrt(20)/.2
    power = stats.nct.cdf(-critical, 19, ncp)+stats.nct.sf(critical, 19, ncp)
    assert power == pytest.approx(.8, abs=1e-10)
    assert .13 < s.planning_mde < .14
    assert s.mde_status == 'model_based_planning'
    assert s.mde_attainable
    assert 'model_based_power_warning' in s.warning_codes


def test_mde_depends_on_planning_assumptions_not_observed_effect():
    one = call(np.linspace(-.3, .5, 20), planning_sd=.2).summary().iloc[0]
    shifted = call(np.linspace(.1, .9, 20), planning_sd=.2).summary().iloc[0]
    larger_sd = call(np.linspace(-.3, .5, 20), planning_sd=.4).summary().iloc[0]
    more_patients = call(np.linspace(-.3, .5, 40), planning_sd=.2).summary().iloc[0]
    higher_power = call(np.linspace(-.3, .5, 20), planning_sd=.2, planning_power=.9).summary().iloc[0]
    assert one.planning_mde == shifted.planning_mde
    assert larger_sd.planning_mde == pytest.approx(2*one.planning_mde)
    assert more_patients.planning_mde < one.planning_mde < higher_power.planning_mde


def test_mde_not_clipped_to_possible_correlation_scale():
    s = call([-.6, -.2, .1, .2, .4, .6], planning_sd=1.).summary().iloc[0]
    assert s.planning_mde > 1
    assert not s.mde_attainable
    assert 'mde_exceeds_correlation_range' in s.warning_codes


@pytest.mark.parametrize('j,alpha,power', [(2, .05, .8), (3, .05, .8),
                                           (4, .05, .8), (6, .05, .99),
                                           (4, .9, .95)])
def test_small_j_mde_matches_independent_normal_chi_integral(j, alpha, power):
    s = call(np.linspace(-.3, .5, j), min_patients=2, planning_sd=.2,
             alpha=alpha, planning_power=power).summary().iloc[0]
    assert s.mde_status == 'model_based_planning'
    df = j - 1
    critical = stats.t.isf(alpha/2, df)
    delta = s.planning_mde*np.sqrt(j)/.2
    def integrand(chi):
        boundary = critical*chi/np.sqrt(df)
        return (special.ndtr(delta-boundary)+special.ndtr(-delta-boundary))*stats.chi.pdf(chi, df)
    attained, error = integrate.quad(integrand, 0, np.inf, epsabs=1e-10, epsrel=1e-10)
    assert error < 1e-9
    assert attained == pytest.approx(power, abs=1e-9)


def test_planning_source_is_recorded_as_user_declared_not_verified():
    values = np.linspace(-.3, .5, 8)
    unspecified = call(values, planning_sd=.2)
    specified = call(values, planning_sd=.2, planning_sd_source=' independent pilot cohort ')
    assert 'planning_sd_source_unspecified' in unspecified.summary().iloc[0].warning_codes
    assert 'planning_sd_source_unspecified' not in specified.summary().iloc[0].warning_codes
    assert specified.metadata['planning_sd_source'] == 'independent pilot cohort'
    assert specified.metadata['planning_sd_source_verified'] is False
    assert specified.summary().iloc[0].planning_mde == unspecified.summary().iloc[0].planning_mde


@pytest.mark.parametrize('source', [True, 2, '', '   ', np.nan])
def test_planning_source_must_be_nonempty_text(source):
    with pytest.raises(ValueError, match='planning_sd_source'):
        call(np.linspace(-.3, .5, 8), planning_sd=.2, planning_sd_source=source)


def test_planning_source_requires_sd():
    with pytest.raises(ValueError, match='planning_sd'):
        call(np.linspace(-.3, .5, 8), planning_sd_source='External pilot')


def test_external_planning_is_distinct_from_degenerate_observed_interval():
    out = call([.2] * 8, planning_sd=.2)
    s = out.summary().iloc[0]
    assert s.status == 'withheld_technical'
    assert not s.inference_available
    assert pd.isna(s.ci_low)
    assert s.mde_status == 'model_based_planning'
    assert np.isfinite(s.planning_mde)


@pytest.mark.parametrize('aggregate', ['re_hksj', 'eb'])
def test_mde_never_reuses_equal_scale_for_other_targets(aggregate):
    res = call(np.linspace(-.3, .5, 8), aggregate=aggregate, planning_sd=.2)
    s = res.summary().iloc[0]
    assert pd.isna(s.planning_mde)
    assert s.mde_status == 'not_supported'
    assert 'mde_not_supported' in s.warning_codes
    assert 'weighted_aggregation_assumptions' in s.warning_codes


def test_no_planning_preserves_valid_high_alpha():
    s = call(np.linspace(-.3, .5, 8), alpha=.9).summary().iloc[0]
    assert s.status == 'reported'
    assert s.confidence_level == pytest.approx(.1)
    assert s.mde_status == 'not_requested'


def test_requested_planning_power_must_exceed_alpha():
    with pytest.raises(ValueError, match='planning_power'):
        call(np.linspace(-.3, .5, 8), alpha=.9, planning_sd=.2)


def test_eb_below_floor_has_no_individual_interval_bypass():
    res = call([-.5, -.2, .1, .3, .6], aggregate='eb', planning_sd=.2)
    s = res.summary().iloc[0]
    assert s.reporting_status == 'descriptive_only'
    assert s.mde_status == 'withheld_policy'
    assert res.patients[['eb_estimate', 'eb_ci_low', 'eb_ci_high']].isna().all().all()


@pytest.mark.parametrize('keyword,value', [
    ('max_ci_width', True), ('max_ci_width', np.nan), ('max_ci_width', np.inf),
    ('max_ci_width', 0), ('max_ci_width', -.2), ('max_ci_width', '0.2'),
    ('planning_sd', False), ('planning_sd', np.nan), ('planning_sd', np.inf),
    ('planning_sd', 0), ('planning_sd', 1.01),
    ('minimum_relevant_effect', True), ('minimum_relevant_effect', np.nan),
    ('minimum_relevant_effect', np.inf), ('minimum_relevant_effect', -.1),
    ('minimum_relevant_effect', 0), ('minimum_relevant_effect', 1.01),
    ('planning_power', True), ('planning_power', np.nan), ('planning_power', np.inf),
    ('planning_power', 1), ('planning_power', .04), ('planning_power', 0)])
def test_invalid_reporting_parameters_fail(keyword, value):
    kwargs = {keyword: value}
    if keyword == 'minimum_relevant_effect':
        kwargs['planning_sd'] = .2
    if keyword == 'planning_power':
        kwargs['planning_sd'] = .2
    with pytest.raises(ValueError):
        call(np.linspace(-.3, .5, 8), **kwargs)


def test_relevant_effect_requires_independent_sd():
    with pytest.raises(ValueError, match='planning_sd'):
        call(np.linspace(-.3, .5, 8), minimum_relevant_effect=.1)


def test_plots_do_not_draw_withheld_cohort_inference_and_keep_cautions():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out = call([-.3, .1, .4, .5, .6])
    ax = out.plot_forest(title='Custom title')
    assert len(ax.get_yticklabels()) == 5
    assert 'withheld' in ax.get_title().lower()
    assert all(line.get_marker() != 'D' for line in ax.lines)
    plt.close(ax.figure)
    out = call(np.linspace(-.3, .5, 8))
    ax = out.plot_forest(title='Custom title')
    assert 'caution' in ax.get_title().lower()
    plt.close(ax.figure)
