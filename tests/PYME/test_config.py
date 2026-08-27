import pytest
from unittest.mock import patch, MagicMock


def _make_ep(name, value):
    ep = MagicMock()
    ep.name = name
    ep.value = value
    return ep


def _fresh_parse():
    """Re-run _parse_plugin_config against a clean plugins dict."""
    import PYME.config as cfg
    cfg.plugins = {app: set() for app in ['visgui', 'dsviewer', 'recipes', 'fit_factories']}
    cfg._parse_plugin_config()
    return cfg


class TestEntryPointDiscovery:
    def test_recipe_entry_point_added(self):
        eps = {'pyme.plugins.recipes': [_make_ep('myplugin', 'mypackage.recipe_modules')]}

        def fake_entry_points(group):
            return eps.get(group, [])

        with patch('PYME.config.entry_points', side_effect=fake_entry_points, create=True):
            cfg = _fresh_parse()

        assert 'mypackage.recipe_modules' in cfg.get_plugins('recipes')

    def test_visgui_entry_point_added(self):
        eps = {'pyme.plugins.visgui': [_make_ep('myplugin', 'mypackage.visgui_modules')]}

        def fake_entry_points(group):
            return eps.get(group, [])

        with patch('PYME.config.entry_points', side_effect=fake_entry_points, create=True):
            cfg = _fresh_parse()

        assert 'mypackage.visgui_modules' in cfg.get_plugins('visgui')

    def test_colon_attr_suffix_stripped(self):
        """Only the module path before ':' should be used."""
        eps = {'pyme.plugins.recipes': [_make_ep('myplugin', 'mypackage.recipe_modules:some_attr')]}

        def fake_entry_points(group):
            return eps.get(group, [])

        with patch('PYME.config.entry_points', side_effect=fake_entry_points, create=True):
            cfg = _fresh_parse()

        assert 'mypackage.recipe_modules' in cfg.get_plugins('recipes')
        assert 'mypackage.recipe_modules:some_attr' not in cfg.get_plugins('recipes')

    def test_multiple_groups_and_entries(self):
        eps = {
            'pyme.plugins.recipes': [
                _make_ep('a', 'pkg.recipes_a'),
                _make_ep('b', 'pkg.recipes_b'),
            ],
            'pyme.plugins.fit_factories': [_make_ep('c', 'pkg.fitters')],
        }

        def fake_entry_points(group):
            return eps.get(group, [])

        with patch('PYME.config.entry_points', side_effect=fake_entry_points, create=True):
            cfg = _fresh_parse()

        assert {'pkg.recipes_a', 'pkg.recipes_b'} <= cfg.get_plugins('recipes')
        assert 'pkg.fitters' in cfg.get_plugins('fit_factories')

    def test_bad_entry_point_does_not_crash(self):
        """A malformed entry point should log a warning but not abort discovery."""
        bad = MagicMock()
        bad.name = 'bad'
        type(bad).value = property(lambda self: (_ for _ in ()).throw(RuntimeError('broken')))
        good = _make_ep('good', 'pkg.good_module')

        eps = {'pyme.plugins.recipes': [bad, good]}

        def fake_entry_points(group):
            return eps.get(group, [])

        with patch('PYME.config.entry_points', side_effect=fake_entry_points, create=True):
            cfg = _fresh_parse()

        assert 'pkg.good_module' in cfg.get_plugins('recipes')

    def test_no_entry_points_leaves_dict_intact(self):
        """With no entry points declared, existing empty sets remain."""
        def fake_entry_points(group):
            return []

        with patch('PYME.config.entry_points', side_effect=fake_entry_points, create=True):
            cfg = _fresh_parse()

        for app in ['visgui', 'dsviewer', 'recipes', 'fit_factories']:
            assert isinstance(cfg.get_plugins(app), set)
