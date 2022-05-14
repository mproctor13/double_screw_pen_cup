#!/bin/bash


OSC=openscad-nightly
PARTS="vase"
SETS="narrow_7_7 wide_4_7 wide_4_4"
NUTS="Left Right Bottom Screw Top"

for SET in $SETS; do 
  echo "Output $SET image"
  $OSC -p double_screw_pen_cup.json -P ${SET} -D "what_to_output=\"model\"" \
    --autocenter --viewall -o img/${SET}_model.png double_screw_pen_cup.scad
  for part in $PARTS; do
    echo "Output ${SET}_${part}.stl"
    $OSC -p double_screw_pen_cup.json -P ${SET} -D "what_to_output=\"${part}\"" \
      --autocenter --viewall -o stl/${SET}_${part}.stl double_screw_pen_cup.scad
  done
  for part in $NUTS; do
    $OSC -p double_screw_pen_cup.json -P ${SET} -D "what_to_output=\"${part}\"" -D "nut_type=\"round\"" \
      -D "nut_radius=40" --autocenter --viewall -o stl/${SET}_${part}_hex.stl double_screw_pen_cup.scad
  done
  for part in $NUTS; do
    $OSC -p double_screw_pen_cup.json -P ${SET} -D "what_to_output=\"${part}\"" -D "nut_type=\"hex\"" \
       -D "nut_radius=44" --autocenter --viewall -o stl/${SET}_${part}_hex.stl double_screw_pen_cup.scad
  done
done

