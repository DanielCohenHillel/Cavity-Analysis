#!/usr/bin/env python
# -*- coding: utf-8 -*-
import gdspy
import numpy as np
import designtools as dt

lib = gdspy.GdsLibrary()      # GDSII file (library)
cell = lib.new_cell('wafer')  # main cell

# =================================👉 Parameters 👈===============================
# Wafer
wafer_rad  = 25400
wafer_del  = np.pi/10  # The angle to begin the wafer (straight cut on the side)
flip_wafer = True      # Flip wafer along the y axis

# Cut details (cut = single saphire chip)
cut_x        = -12900  # x-position of left edge
cut_y_first  = -15250  # y-position of bottom edge of the lowest chip
cut_length   = 25200   # Length of the chip (horizontal)
cut_width    = 1900    # Width of the chip  (vertical)

cut_gap = 300  # Gap between the cuts
cut_num = 14   # Amount of cuts

# RO line details
RO_rel_x   = 12595-1750  # x relative to the bottom left corner of the cut
RO_rel_y   = 0           # y relative to the middle of the cut (should be 0)
RO_length  = 8880
RO_width   = 150

# LC resonator details (zig-zag line)
LC_dist    = 1000   # Distance between LC-resonator and RO-line
LC_width   = 150
LC_height  = 1500  # distance between the zig and the zag (horizontal lines) 👌
LC_rel_y   = 0     # Should be zero
LC_rounds  = 10    # Number of zig-zags
LC_gap     = 300   # Gap between the vertical lines
roundness  = 80    # coil roundness (0~150)

# Layer information (doesn't matter too much)
_conductor = {"layer": 0, "datatype": 0}
_cut       = {"layer": 2, "datatype": 2}
_wafer     = {"layer": 2, "datatype": 3}

# Marking parameters
saww_width = 300  # Width of the saw that cuts the wafer
mark_dist  = 450  # Distance between the markings and the cuts

# Clamp = thing that holds the chip
clamp_length = 3100  # Amount of the chip not inside the cavity

txt_size = 200
txt      = "merd"

# ================================👉 Draw 👈=======================================
# Loop over all the chips
for i in range(cut_num):
    cut_y = cut_y_first + (cut_width+cut_gap)*i  # y of bottom edge of the i-th chip
    # RO position (bottom left corner)
    RO_x  = cut_x + RO_rel_x
    RO_y  = cut_y + cut_width/2 - RO_width/2 + RO_rel_y

    # LC position (bottom left corner)
    LC_x = RO_x  + RO_length + LC_dist
    LC_y = cut_y + cut_width/2 - LC_height/2 + LC_rel_y

    # Create chip rectangle
    cut = gdspy.Rectangle((cut_x, cut_y), (cut_x+cut_length, cut_y + cut_width), **_cut)
    # Create RO rectangle
    RO = gdspy.Rectangle((RO_x, RO_y), (RO_x+RO_length, RO_y+RO_width), **_conductor)

    cell.add(cut)
    cell.add(RO)

    # ======= Draw inductor coil "transmon wannabe"🧲 =======
    path = gdspy.Path(LC_width, (LC_width/2+LC_x, LC_y))
    path.segment(LC_height-LC_width-roundness, "+y")
    for i in range(LC_rounds-1):
        path.turn(LC_width/2+roundness, ["r", "l"][i%2])
        path.segment(LC_gap-roundness*2)
        path.turn(LC_width/2+roundness, ["r", "l"][i%2])
        path.segment(LC_height-2*LC_width+((i+2)//LC_rounds)*(LC_width+roundness) - roundness*2)
    cell.add(path)

    # ========= Rectangular coil ===========
    # for j in range(LC_rounds):
    #     # Vertical lines
    #     rect_x   = LC_x + j*(LC_gap+LC_width)
    #     rect_y   = LC_y
    #     vertical = gdspy.Rectangle((rect_x, rect_y), (rect_x+LC_width, rect_y+LC_height), **_conductor)
    #     cell.add(vertical)

    #     # Bridges between lines
    #     if j==LC_rounds-1:
    #         break;
    #     bridges = gdspy.Rectangle((rect_x+LC_width, rect_y+(LC_height-LC_width)*(1-j%2)),
    #                         (rect_x+LC_width+LC_gap, rect_y+(LC_height-LC_width)*(1-j%2) + LC_width), **_conductor)
    #     cell.add(bridges)
    
    # Dice marks around each cut for the saww
    marks = [
        gdspy.Polygon(dt.dice_mark(cut_x-mark_dist, cut_y + cut_width)),
        gdspy.Polygon(dt.dice_mark(cut_x-mark_dist, cut_y, flip_y=True)),
        gdspy.Polygon(dt.dice_mark(cut_x+cut_length+mark_dist, cut_y + cut_width, flip_x=True)),
        gdspy.Polygon(dt.dice_mark(cut_x+cut_length+mark_dist, cut_y, flip_x=True, flip_y=True))
        ]
    cell.add(marks)

    # Fun text :)
    cell.add(gdspy.Text(f"{txt}/{i+1}", txt_size, angle=np.pi/2, position=(cut_x+500,cut_y+cut_width/2-txt_size*2), layer=3))

    # Triangle marking for clamp insertion depth
    mark_amount = 7
    for i in range (mark_amount):
        gap  = 250
        size = 200
        dx   =  size - (size/mark_amount)*abs(i-mark_amount//2)
        dy   = -size + (size/mark_amount)*abs(i-mark_amount//2)
        x    = cut_x + clamp_length+200*(mark_amount//2-i)
        y    = cut_y + cut_width-dy - gap
        mark = [
            gdspy.Polygon(dt.tri(x,y,dx,dy) , layer=3),  # Top marks
            gdspy.Polygon(dt.tri(x,y-cut_width+2*dy+2*gap,dx,-dy) , layer=3)  # Bottom marks
            ]
        cell.add(mark)


# Ruler
unit_marks = [gdspy.Rectangle([cut_x+600, cut_y_first+100*i-200], [cut_x+700, cut_y_first+15+100*i-200], layer=3) for i in range(310)]
cell.add(unit_marks)

cut_y_last = cut_y_first+cut_num*(cut_width+cut_gap)-cut_gap  # y of top chip
# Marks in the corner of the entire thing
marks = [
    # Bottom left
    gdspy.Polygon(dt.dice_mark(cut_x, cut_y_first-mark_dist, angle=-np.pi/2)),
    gdspy.Polygon(dt.dice_mark(cut_x-saww_width, cut_y_first-mark_dist, angle=-np.pi/2, flip_y=True)),
    # Bottom right
    gdspy.Polygon(dt.dice_mark(cut_x+cut_length+saww_width, cut_y_first-mark_dist, angle=-np.pi/2)),
    gdspy.Polygon(dt.dice_mark(cut_x+cut_length, cut_y_first-mark_dist, angle=-np.pi/2, flip_y=True)),
    # Top left
    gdspy.Polygon(dt.dice_mark(cut_x, cut_y_last + mark_dist, angle=np.pi/2, flip_y=True)),
    gdspy.Polygon(dt.dice_mark(cut_x-saww_width, cut_y_last + mark_dist, angle=np.pi/2)),
    # Top left
    gdspy.Polygon(dt.dice_mark(cut_x+cut_length+saww_width, cut_y_last + mark_dist, angle=np.pi/2, flip_y=True)),
    gdspy.Polygon(dt.dice_mark(cut_x+cut_length, cut_y_last + mark_dist, angle=np.pi/2)),

    # Extra mark around each corner
    gdspy.Polygon(dt.dice_mark(cut_x-mark_dist, cut_y_first - saww_width)),
    gdspy.Polygon(dt.dice_mark(cut_x+cut_length+mark_dist, cut_y_first - saww_width, flip_x=True)),
    gdspy.Polygon(dt.dice_mark(cut_x-mark_dist, cut_y_last + saww_width, flip_y=True)),
    gdspy.Polygon(dt.dice_mark(cut_x+cut_length+mark_dist, cut_y_last + saww_width, flip_x=True, flip_y=True)),
]
cell.add(marks)

# === Draw Wafer ===
wafer_pnts = dt.gen_wafer(wafer_rad, wafer_del, flip_wafer)
cell.add(gdspy.Polygon(wafer_pnts, **_wafer))

# Save to file
lib.write_gds('chips.gds')
# Save as SVG as well
cell.write_svg('chips.svg')

tot_len   = (LC_height+LC_gap)*LC_rounds - LC_gap  # Total length of coil
tot_width = (LC_width +LC_gap)*LC_rounds - LC_gap  # Total width of coil
print('''\n
           Parameters
 ===========================  
| Coil                      |
|   * Coil length: {:.3f}mm |
|   * Coil width:  {:.4f}mm |
|   * Coil height: {:.4f}mm |
| Wafer                     |
|   * Wafer radius: {:.3f}mm|
|   * Cut angle:    {:.4f}  |
| Chips                     |
|   * Chip amount: {:.3f}   |
|   * Chip length: {:.3f}mm |
|   * Chip width:  {:.4f}mm |
| RO line                   |
|   * RO width:   {:.4f}mm  |
|   * RO length:  {:.4f}mm  |
|   * RO-LC dist: {:.4f}mm  |
 ===========================  \n
'''.format(tot_len/1e3, tot_width/1e3, LC_height/1e3, wafer_rad/1e3, wafer_del, cut_num, cut_length/1e3, cut_width/1e3, RO_width/1e3, RO_length/1e3, LC_dist/1e3))

# Display using internal viewer.
gdspy.LayoutViewer()
