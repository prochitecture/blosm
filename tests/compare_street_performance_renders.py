"""Compare full-size controlled renders; optionally create a documentation sheet.

Uses host Python, NumPy and Pillow. Resizing applies only to the presentation
sheet, never to comparison pixels. Generated reports stay in artifacts; an
illustration explicitly included in documentation belongs in the repository.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont


def illustration(reference, candidate, path, identical):
    width, margin, gap = 1240, 24, 16
    tile_width = (width - 2 * margin - gap) // 2
    tile_height = tile_width * 2 // 3
    canvas = Image.new('RGB', (width, 180 + 2 * (tile_height + 52) + 65), '#eef2f5')
    draw = ImageDraw.Draw(canvas)
    font_path = Path('C:/Windows/Fonts/segoeui.ttf')
    font = lambda size: ImageFont.truetype(str(font_path), size) if font_path.exists() else ImageFont.load_default()
    draw.text((margin, 17), 'Street performance: preserved appearance', font=font(32), fill='#182733')
    text = ('Eight paired renders: no changed pixels at full resolution.' if identical else
            'Render differences detected; consult the full-resolution comparison report.')
    draw.text((margin, 65), text, font=font(20), fill='#17614b' if identical else '#a33131')
    for column, label in enumerate(('Reference', 'Optimized')):
        draw.text((margin + column * (tile_width + gap), 116), label, font=font(25), fill='#182733')
    for row, (name, caption) in enumerate((('junction_relief.png', 'Curved surface · exaggeration 10'),
                                          ('junction_wire.png', 'Direct quads · rounded boundaries'))):
        y = 162 + row * (tile_height + 52)
        for column, directory in enumerate((reference, candidate)):
            with Image.open(directory / name) as image:
                canvas.paste(image.convert('RGB').resize((tile_width, tile_height), Image.Resampling.LANCZOS),
                             (margin + column * (tile_width + gap), y))
        draw.text((margin, y + tile_height + 8), caption, font=font(20), fill='#415462')
    draw.text((margin, canvas.height - 45), 'Identical cameras, lighting and shaders · CPU Cycles · 32 samples · no denoising',
              font=font(18), fill='#415462')
    path.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--illustration', type=Path)
    args = parser.parse_args()
    a, b = [json.loads((p / 'renders.json').read_text(encoding='utf8'))
            for p in (args.reference, args.candidate)]
    for key in ('samples', 'seed', 'mesh_count', 'shots'):
        assert a[key] == b[key], ('Render setup differs', key)
    assert len(a['shots']) == 8, 'Require the complete controlled-view set'
    records = []
    for shot in a['shots']:
        name = shot['image']
        with Image.open(args.reference / name) as first, Image.open(args.candidate / name) as second:
            assert first.size == second.size and first.mode == second.mode
            delta = np.abs(np.asarray(first).astype(np.int16) - np.asarray(second).astype(np.int16))
        records.append(dict(image=name, maximum_channel_difference=int(delta.max()),
                            changed_pixels=int(np.any(delta, axis=-1).sum()),
                            mean_absolute_difference=float(delta.mean())))
    passed = all(row['changed_pixels'] == 0 for row in records)
    report = dict(passed=passed, shots=records, settings_identical=True,
                  full_resolution=True, samples=a['samples'], seed=a['seed'])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding='utf8')
    if args.illustration:
        illustration(args.reference, args.candidate, args.illustration, passed)
    print(json.dumps(report), flush=True)
    assert passed, 'Controlled render pixels differ'


if __name__ == '__main__':
    main()
