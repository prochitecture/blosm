"""Register/unregister street runtime in Blender's actual restricted context."""
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import bpy
from _bpy_restrict_state import RestrictBlend

sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from blosm.way.gn.editable_input import model,runtime


class StartupTests(unittest.TestCase):
    def setUp(self):
        model.register()

    def tearDown(self):
        runtime.unregister()
        model.unregister()

    def test_restricted_registration_defers_access_until_data_exists(self):
        with patch.object(runtime,'_index_legacy_sources',wraps=runtime._index_legacy_sources) as index:
            with RestrictBlend():
                runtime.register();runtime.register()
                self.assertFalse(hasattr(bpy.data,'objects'))
                self.assertEqual(index.call_count,0)
                self.assertTrue(bpy.app.timers.is_registered(runtime._initialize_sources))
                self.assertEqual(runtime._initialize_sources(),.05)
            self.assertIsNone(runtime._initialize_sources())
            self.assertEqual(index.call_count,1)
        self.assertEqual(bpy.app.handlers.load_post.count(runtime._loaded),1)
        self.assertEqual(bpy.app.handlers.depsgraph_update_post.count(runtime._changed),1)
        runtime.unregister()
        self.assertFalse(bpy.app.timers.is_registered(runtime._initialize_sources))
        self.assertFalse(bpy.app.timers.is_registered(runtime._tick))

    def test_disable_before_deferred_initialization_leaves_no_callback(self):
        with RestrictBlend():
            runtime.register();runtime.unregister()
        self.assertFalse(bpy.app.timers.is_registered(runtime._initialize_sources))
        self.assertIsNone(runtime._initialize_sources())

    def test_normal_registration_indexes_immediately_and_load_rebuilds_once(self):
        with patch.object(runtime,'_index_legacy_sources',wraps=runtime._index_legacy_sources) as index:
            runtime.register();self.assertEqual(index.call_count,1)
            self.assertFalse(bpy.app.timers.is_registered(runtime._initialize_sources))
            runtime._loaded(None);self.assertEqual(index.call_count,2)

    def test_complete_addon_registration_in_restricted_context(self):
        import blosm
        try:
            with RestrictBlend():blosm.register()
            self.assertTrue(bpy.app.timers.is_registered(runtime._initialize_sources))
            self.assertIsNone(runtime._initialize_sources())
        finally:
            blosm.unregister()


if __name__=='__main__':
    result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(StartupTests))
    if not result.wasSuccessful():raise AssertionError('Street runtime startup regression failed')
