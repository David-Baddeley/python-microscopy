.. _plugins:

Plugin Registration and Discovery
**********************************

PYME supports four plugin groups:

.. csv-table::
    :header: Group, Used by, Plug() argument type
    :widths: auto

    ``visgui``,       PYMEVisualise (VisGUI),  :class:`~PYME.LMVis.VisGUI.VisGUIFrame`
    ``dsviewer``,     PYMEImage / dh5view,     :class:`~PYME.DSView.dsviewer.DSViewFrame`
    ``recipes``,      Recipes / all apps,      *(imported for side-effects; no Plug)*
    ``fit_factories``,Localisation analysis,   *(imported for side-effects; no Plug)*

Plugins in all groups are discovered from the sources below and merged — in every case the module is *not* imported
during discovery, only at the point the relevant application component initialises.


importlib.metadata entry points (recommended)
=============================================

Declare entry points in your package's ``pyproject.toml``. PYME discovers them automatically on startup with no
file copying. This works identically for pip, conda, and editable installs.

.. code-block:: toml

    [project.entry-points."pyme.plugins.recipes"]
    myplugin        = "mypackage.recipe_modules"
    myplugin_extra  = "mypackage.extra_modules"

    [project.entry-points."pyme.plugins.visgui"]
    myplugin = "mypackage.visgui_modules"

    [project.entry-points."pyme.plugins.dsviewer"]
    myplugin = "mypackage.dsviewer_modules"

    [project.entry-points."pyme.plugins.fit_factories"]
    myplugin = "mypackage.fit_factories"

The entry point *name* (left-hand side) is arbitrary and ignored by PYME; the *value* is the fully qualified module
path. An optional ``:attr`` suffix (e.g. ``mypackage.module:some_func``) is accepted but only the module path is used.


YAML config files (required for report plugins; otherwise legacy)
=================================================================

YAML config files are the only current mechanism for registering ``reports`` plugins (templates and filters).
For the four main plugin groups, prefer entry points instead.

Drop a ``<plugin-name>.yaml`` file into any PYME config directory
(``~/.PYME/plugins/``, ``/etc/PYME/plugins/``, or ``<sys.prefix>/etc/PYME/plugins/``).

.. code-block:: yaml

    # visgui/dsviewer/recipes/fit_factories sections work but entry points are preferred.

    reports:
        templates: mypackage.report_templates
        filters:
            mypackage.report_filters:
                - myfilter1

See :func:`PYME.config.get_plugins` and :mod:`PYME.config` for config directory locations.


Legacy .txt files
=================

Individual ``<name>.txt`` files placed in ``plugins/visgui/``, ``plugins/dsviewer/``, or ``plugins/recipes/``
subdirectories of any config directory are still supported.  Each line of the file should be a fully qualified module
path.  New plugins should use one of the mechanisms above instead.


Writing plugins
===============

``visgui`` and ``dsviewer`` plugins must implement a top-level ``Plug(parent)`` function.  ``recipes`` and
``fit_factories`` modules self-register during import via decorator or class-level code — no ``Plug`` is needed.

See :ref:`extendingdsviewer`, :ref:`extendingvisgui`, and :ref:`writingrecipemodules` for details.
