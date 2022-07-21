# Strained-Moire-Visualization
This script library is for the visualization of Moire lattices under a universal 2D strain as described in  URL. The scripts were written in the mathematica 11.3 version. For any use of the figures created by these scripts cite LINK TO THE MANUSCRIPT.

# Basics of the script
The script is used to make images of a Moire under strain for various underlying monolayer lattices. The script generates the unit cells of the monolayers and then overlays them such that a Moire image is generated. The images are continuously adjustable for several different twist angles and strain parameters. The moire vectors are calculated according to the mathematical algorithm LINK TO THE MANUSCRIPT. Then the Moire lattice is shown as 4 red parallellograms and overlayed on top of the Moire lattice such that the corners of the parallellograms point to the commensurate points (where the atoms are aligned on top of each other). The Moire image and the parallellograms are calculated separately. The parallellograms follow the commensurate points of the moire image perfectly for all settings, thus proving the correctness of the mathematical algorithm in LINK TO THE MANUSCRIPT.

![image](https://user-images.githubusercontent.com/61543540/177525456-2112e281-a809-464d-9249-3afa660a1177.png)

# Usage of the script
1)  Compile the cell for the Bravais lattice
2)  Compile the cell for the functions to calculate the moire lattice
3)  Compile the cell to visualize the moire lattice
Use the manipulate sliders to adjust the moire lattice to the wanted values of strain and parameters. The "Amount of unit cells shown" slider adjusts the number of unit cells used for the image. For large Moire patterns this setting has to be large which increases computation time.

# Create a new bilayer moire
Use the "Create_new_bilayer" file. To create a new bilayer moire it's only necessary to change the first cell of the script. 
1)  Set the right parameters for the Bravas lattice
```
a := {0.3498, 0.6338} (*Type in lattice constants of bilayers*)
beta := 90 (*Type in angle between the primitive lattice vectors*)
delta := 0.0(*Adjust the lattice mismatch for your bilayer structure*)
```
2)  Create all the atomes in the primitive cell
```
unitWpos1[x_, y_] := {a[[1]]*x, a[[2]]*y}
unitWpos2[x_, y_] := {a[[1]]*(x + 0.5), a[[2]] (y + 0.4)}
unitTepos1[x_, y_] := {a[[1]]*x, a[[2]]*(y + 0.3)}
unitTepos2[x_, y_] := {a[[1]]*x, a[[2]]*(y + 0.65)}
unitTepos3[x_, y_] := {a[[1]]*(x + 0.5), a[[2]]*(y + 0.14)}
unitTepos4[x_, y_] := {a[[1]]*(x + 0.5), a[[2]]*(y + 0.79)}
nextpos1[x_, y_] := {a[[1]]*(x + 1), a[[2]]*(y + 0.3)}
nextpos2[x_, y_] := {a[[1]]*(x + 1), a[[2]]*(y + 0.65)}
nextpos3[x_, y_] := {a[[1]]*(x + 1), a[[2]]*(y + 0.)}
nextpos4[x_, y_] := {a[[1]]*(x + 1), a[[2]]*(y + 1)}
nextpos5[x_, y_] := {a[[1]]*(x + 0), a[[2]]*(y + 1)}
nextpos6[x_, y_] := {a[[1]]*(x + 0), a[[2]]*(y - 0.3)}
```
3)  Create all the bonds between the different atoms in the unit cell and to next unit cells
```
unitCell[x_, y_] := {Gray,
  Line[{unitWpos1[x, y], unitTepos1[x, y]}], 
  Line[{unitWpos1[x, y], nextpos6[x, y]}], 
  Line[{unitWpos1[x, y], unitTepos3[x, y]}], 
  Line[{unitWpos2[x, y], unitTepos3[x, y]}], 
  Line[{unitWpos2[x, y], unitTepos1[x, y]}], 
  Line[{unitWpos2[x, y], unitTepos2[x, y]}], 
  Line[{unitWpos2[x, y], unitTepos4[x, y]}], 
  Line[{unitWpos2[x, y], nextpos1[x, y]}], 
  Line[{unitWpos2[x, y], nextpos2[x, y]}], 
  Line[{unitTepos3[x, y], nextpos3[x, y]}], 
  Line[{unitTepos4[x, y], nextpos4[x, y]}], 
  Line[{unitTepos4[x, y], nextpos5[x, y]}], Blue, 
  Disk[unitWpos1[x, y], a[[1]]*0.1], 
  Disk[unitWpos2[x, y], a[[1]]*0.1], Red, 
  Disk[unitTepos1[x, y], a[[1]]*0.1], 
  Disk[unitTepos2[x, y], a[[1]]*0.1], 
  Disk[unitTepos3[x, y], a[[1]]*0.1], 
  Disk[unitTepos4[x, y], a[[1]]*0.1]}
```
4)  Check that your script generates a proper lattice
```
Graphics[Block[{unitVectA = {1, 0}, 
   unitVectB = {Cos[90 Degree], Sin[90 Degree]}}, 
  Table[unitCell @@ (unitVectA j + unitVectB k), {j, -2, 
    2}, {k, -2 + Ceiling[j], 2 + Ceiling[j]}]], 
 PlotRange -> {{-1*a[[2]], 1*a[[2]]}, {-1*a[[2]], 1*a[[2]]}}]
```
5)  Continue compiling the next cells and visualize the custom made moire lattice
