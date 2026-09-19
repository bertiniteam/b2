# This file is part of Bertini 2.
#
# python/test/docs/test_config_reference.py is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# Bertini 2 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.
# See the GNU General Public License for more details.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license, as well as COPYING.

"""The configuration reference page is generated, so what is under test is the generator.

The page (``python/docs/source/detailed/configuration.rst``) is a few paragraphs of prose plus
one directive; every setting, default and description on it is read out of the library when the
docs are built.  That is what b2#406 asked for -- a hand-maintained table of defaults drifts from
the code silently -- but it moves the failure: instead of the page going stale, the generator can
go quietly wrong, printing nothing, or printing a default nobody can trust.

So these check the generator against the library rather than against a golden page: every config
is found, every field has a default worth printing, and a handful of values are what the code
says they are.  They do not run Sphinx; the docs build does that, and does it with warnings as
errors.
"""

import importlib
import os
import sys

import pytest

_EXT = os.path.join(os.path.dirname(__file__), '..', '..', 'docs', 'source', '_ext')
sys.path.insert(0, os.path.abspath(_EXT))

config_reference = pytest.importorskip('config_reference')


def _all_rows():
    """Every (config class name, field, default, description) the page would print."""
    out = []
    for module_name, _, _ in config_reference._MODULES:
        module = importlib.import_module(module_name)
        for class_name, cls in config_reference._config_classes(module):
            rows, problem = config_reference._rows_for(cls)
            assert problem is None, f'{class_name}: {problem}'
            for field, default, description in rows:
                out.append((class_name, field, default, description))
    return out


def test_every_config_class_is_found():
    found = set()
    for module_name, _, _ in config_reference._MODULES:
        module = importlib.import_module(module_name)
        found.update(name for name, _ in config_reference._config_classes(module))

    # the ones that existed when the page was written; the point of generating is that NEW
    # ones appear without anybody editing a list, so this is a floor and not an inventory
    expected = {
        'AMPConfig', 'FixedPrecisionConfig', 'NewtonConfig', 'SteppingConfig',
        'CauchyConfig', 'EndgameConfig', 'PowerSeriesConfig', 'SecurityConfig',
        'AutoRetrackConfig', 'MidPathConfig', 'PostProcessingConfig', 'RecordsConfig',
        'RegenerationConfig', 'SharpeningConfig', 'TolerancesConfig', 'ZeroDimConfig',
    }
    assert expected <= found, f'missing from the reference: {sorted(expected - found)}'


def test_every_field_has_a_default_worth_printing():
    # A field whose default cannot be read is the failure this page exists to prevent: a
    # reference that says nothing, or worse, prints whatever an uninitialized member held.
    unreadable = [(c, f) for c, f, d, _ in _all_rows() if 'unavailable' in d]
    assert not unreadable, f'no default could be read for: {unreadable}'


def test_every_field_is_described():
    undocumented = [(c, f) for c, f, _, desc in _all_rows() if desc == '*undocumented*']
    assert not undocumented, f'no docstring for: {undocumented}'


def test_the_settings_surface_and_the_page_agree():
    # The page lists what update() accepts, because both read writable_fields.  If they ever
    # disagreed the page would be documenting settings nobody can set, or missing ones they can.
    from bertini.config import writable_fields
    import bertini.nag_algorithm as na

    listed = {f for c, f, _, _ in _all_rows() if c == 'ZeroDimConfig'}
    assert listed == set(writable_fields(na.ZeroDimConfig))


@pytest.mark.parametrize('config_name,field,expected', [
    ('SharpeningConfig', 'sharpendigits', '``0``'),               # 0 = do not sharpen
    ('EndgameConfig', 'sample_factor', '``1/2``'),                # an exact rational, printed exactly
    ('SteppingConfig', 'max_num_steps', '``100000``'),
    ('SecurityConfig', 'level', '``0``'),
    ('RecordsConfig', 'recall', '``RecallPolicy.Completed``'),    # an enum, printed as one writes it
    ('RegenerationConfig', 'slice_final_tolerance', '``1e-11``'),
])
def test_spot_defaults(config_name, field, expected):
    rows = {(c, f): d for c, f, d, _ in _all_rows()}
    assert rows[(config_name, field)] == expected


def test_the_adaptive_bounds_say_they_come_from_the_system():
    # A default-constructed AMPConfig leaves these unset on purpose -- they are derived from the
    # system by amp_config_from(system) -- so the page must say that rather than print a number.
    rows = {(c, f): d for c, f, d, _ in _all_rows()}
    for field in ('degree_bound', 'coefficient_bound', 'jacobian_eval_error_bound'):
        assert rows[('AMPConfig', field)] == '*derived from the system*'
