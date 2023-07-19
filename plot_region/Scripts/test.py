from PIL import Image, ImageDraw, ImageFont
import pandas as pd
import math

# Create some example data
m = pd.DataFrame([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.7, 0.8, 0.9]])
# Add colorbar
L= 3030
W= 1630
margins = {'left':50, 'right':100, 'top':50, 'bottom':50}
# Initiate color pallete
pal = [(0, 85, 63),(0, 132, 98),(0, 180, 132),(0, 227, 167),(19, 255, 193),
       (66, 255, 205),(113, 255, 218),(161, 255, 230),(208, 255, 243),(255, 255, 255)]
palr = list(reversed(pal))
# Initiate image
im = Image.new('RGB', (L, W), (255, 255, 255))
draw = ImageDraw.Draw(im)

colorbar_width = 20
colorbar_height = 200
colorbar_margin = 20
colorbar = Image.new('RGB', (colorbar_width, colorbar_height), (255, 255, 255))
colorbar_draw = ImageDraw.Draw(colorbar)
# Get the minimum and maximum values from m
m_min = m.min().min()
m_max = m.max().max()
# Draw the colorbar
for i in range(colorbar_height):
    v = (colorbar_height - i - 1) / (colorbar_height - 1) * (m_max - m_min) + m_min
    ic = math.floor((v - m_min) / (m_max - m_min) * 9)
    colorbar_draw.line((0, i, colorbar_width, i), fill=palr[ic])
    if i % 50 == 0: # Add ticks and labels every 50 pixels
      tick_label = round(v, 2)
      font = ImageFont.truetype("Ubuntu-R.ttf", 15)
      colorbar_draw.text((colorbar_width + 5, i - 5), str(tick_label), fill=(0, 0, 0), font=font)
      colorbar_draw.line((colorbar_width, i, colorbar_width + 5, i), fill=(0, 0, 0))
# Save the colorbar image to a file
colorbar.save("colorbar.png")
# Paste the colorbar into the image
draw.bitmap((margins['left'] + L + colorbar_margin, margins['top']), colorbar)
