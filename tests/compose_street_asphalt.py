"""Create labeled comparison sheets from the native, unretouched Cycles renders."""
import argparse
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    font = ImageFont.truetype('C:/Windows/Fonts/segoeui.ttf', 20)
    small = ImageFont.truetype('C:/Windows/Fonts/segoeui.ttf', 16)
    for view in ('detail', 'long'):
        width, height = (600, 400) if view == 'detail' else (800, 213)
        canvas = Image.new('RGB', (2 * width + 60, 3 * (height + 65) + 125), '#17212b')
        draw = ImageDraw.Draw(canvas)
        draw.text((20, 12), 'POLY HAVEN ASPHALT — GENERATED CARRIAGEWAYS', font=font, fill='white')
        draw.text((20, 45), 'PBR', font=font, fill='#b5dbef')
        draw.text((40 + width, 45), 'DIFFUSE ONLY', font=font, fill='#b5dbef')
        for row, (wear, source, tile) in enumerate((('NEW', 'Clean Asphalt', '2.10 m'),
                                                   ('WORN', 'Asphalt 01', '2.08 m'),
                                                   ('CRACKED', 'Asphalt 02', '3.00 m'))):
            y = 80 + row * (height + 65)
            for column in range(2):
                index = row * 2 + column
                x = 20 + column * (width + 20)
                draw.text((x, y), f'{wear} · {source} · tile {tile}', font=small, fill='white')
                image = Image.open(args.output / f'{index:02d}_{view}.png').convert('RGB')
                image = image.resize((width, height), Image.Resampling.LANCZOS)
                canvas.paste(image, (x, y + 30))
        draw.text((20, canvas.height - 30),
                  'Same lighting and exposure · 160 m streets · original 2K maps · no displacement or texture retouching',
                  font=small, fill='#b5c4ce')
        canvas.save(args.output / f'asphalt_{view}_comparison.png')
    print('ASPHALT COMPARISON SHEETS READY', args.output)


if __name__ == '__main__':
    main()
