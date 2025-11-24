# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2022 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2022 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2022 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/main/LICENSE

import os
import tempfile
import unittest
from pathlib import Path

from PIL import Image


# Import the function directly to avoid mpas_analysis.test dependencies
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))
from mpas_analysis.shared.html.image_xml import _generate_thumbnails


class TestThumbnailGeneration(unittest.TestCase):
    """Test thumbnail generation with Pillow"""

    def test_generate_thumbnails_horizontal(self):
        """Test thumbnail generation for horizontal images"""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create a horizontal test image
            test_image = Image.new('RGB', (1200, 600), color='blue')
            image_filename = 'test_horizontal.png'
            image_path = Path(tmpdir) / image_filename
            test_image.save(image_path)

            # Generate thumbnails
            imageSize, thumbnailSize, orientation = _generate_thumbnails(
                image_filename, tmpdir
            )

            # Verify results
            self.assertEqual(imageSize, (1200, 600))
            self.assertEqual(orientation, 'horiz')
            self.assertEqual(thumbnailSize[1], 120)  # height should be 120
            
            # Check thumbnail files exist
            thumbnail_dir = Path(tmpdir) / 'thumbnails'
            self.assertTrue(thumbnail_dir.exists())
            self.assertTrue((thumbnail_dir / 'test_horizontal.jpg').exists())
            self.assertTrue((thumbnail_dir / 'fixed_test_horizontal.jpg').exists())

    def test_generate_thumbnails_vertical(self):
        """Test thumbnail generation for vertical images"""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create a vertical test image
            test_image = Image.new('RGB', (400, 800), color='green')
            image_filename = 'test_vertical.png'
            image_path = Path(tmpdir) / image_filename
            test_image.save(image_path)

            # Generate thumbnails
            imageSize, thumbnailSize, orientation = _generate_thumbnails(
                image_filename, tmpdir
            )

            # Verify results
            self.assertEqual(imageSize, (400, 800))
            self.assertEqual(orientation, 'vert')
            self.assertEqual(thumbnailSize[1], 320)  # height should be 320
            
            # Check thumbnail files exist
            thumbnail_dir = Path(tmpdir) / 'thumbnails'
            self.assertTrue(thumbnail_dir.exists())
            self.assertTrue((thumbnail_dir / 'test_vertical.jpg').exists())
            self.assertTrue((thumbnail_dir / 'fixed_test_vertical.jpg').exists())

    def test_image_lanczos_constant(self):
        """Test that Image.LANCZOS constant is available"""
        # This test ensures that Image.LANCZOS is available across
        # Pillow versions 10.x, 11.x, and 12.x
        self.assertTrue(hasattr(Image, 'LANCZOS'))
        self.assertIsNotNone(Image.LANCZOS)
        
        # Test that resize works with LANCZOS
        test_image = Image.new('RGB', (100, 100), color='red')
        resized = test_image.resize((50, 50), Image.LANCZOS)
        self.assertEqual(resized.size, (50, 50))
