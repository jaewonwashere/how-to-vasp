# Example Workflow

These instructions are more or less my workflow for running convergence tests. My POSCAR files come from OCP (Open Catalyst Project) calculations, and I have them all in a directory called poscarfiles. My POSCAR files are named like "POSCAR_mp-126_N", which corresponds to Nitrogen placed onto Pt (mp-126). 

## Prepare the Working Directory

Create a project directory:
```
cd ~
mkdir convergence_(material name)
cd convergence_(material name)
mkdir ktest encut
```

To continue my Pt / N example, I would name this directory "convergence_Pt"

## Get the POSCAR File

If the POSCAR file is named like POSCAR_mp-126_N:
```
cp ~/poscarfiles/POSCAR_mp-126_N ./
mv POSCAR_mp-126_N POSCAR
```

For example, this takes POSCAR_mp-126_N from my poscarfiles directory, moves it to the convergence_Pt directory, and renames it to POSCAR.  

## Set up POTCAR

Once the POSCAR file is renamed and in the convergence_(material name) directory, we can set up an automation script to generate the full POTCAR file.

First, create the following automation script in your home directory. 
```
cd ~
nano generate_POTCAR
```
Then, paste the following code into the script. This automates the process of retrieving and combining the individual POTCAR_(atom) files. 

```
import os
# Input filenames
poscar_path = "POSCAR"
potcar_dir = os.path.expanduser("~/vasp_support/potpaw_PBE")
output_potcar = "POTCAR"
def get_element_order_from_poscar(poscar_file):
    with open(poscar_file, "r") as f:
        lines = f.readlines()
    # Line 6: element symbols
    elements = lines[5].split()
    return elements
def concatenate_potcars(elements, potcar_root, output_file):
    with open(output_file, "wb") as outfile:
        for el in elements:
            el_path = os.path.join(potcar_root, el, "POTCAR")
            if not os.path.exists(el_path):
                raise FileNotFoundError(f"Missing POTCAR for element: {el} at {el_path}")
            with open(el_path, "rb") as potcar_part:
                outfile.write(potcar_part.read())
    print(f"[✓] POTCAR created with elements: {' '.join(elements)}")
if __name__ == "__main__":
    elements = get_element_order_from_poscar(poscar_path)
    concatenate_potcars(elements, potcar_dir, output_potcar)
```

Finally, copy the script from your home directory to your convergence_(material name) directory and run the script to generate the combined POTCAR file. In future runs you will only need to do this step, as the generate_POTCAR will already be in your home directory. 

```
cp ~/generate_POTCAR ./
python generate_POTCAR
```

Ensure the ordering matches the atomic order in the POSCAR by using the following commmand. In our example case, we would want to make sure that the output first displays Pt, then **after** Pt displays N. 

```
grep VRHFIN POTCAR
```


## Import Input Templates

Import all the following files and scripts into your convergence_(material name) directory

```
cp ~/vasp_inputs/templates/INCAR .
cp ~/vasp_inputs/templates/job .
cp ~/vasp_inputs/templates/runKPOINTS.py .
cp ~/vasp_inputs/templates/extractDataKPOINTS.py .
cp ~/vasp_inputs/templates/runENCUT.py .
cp ~/vasp_inputs/templates/extractDataENCUT.py .
```

## Edit INCAR for Convergence Testing

``` 
SYSTEM = [Adsorbate] on [Surface](111)
ENCUT = 450
GGA = PE
# NELM = 200

# scf
# PREC = Accurate
EDIFF = 1e-6 # energy convergence criterion
ALGO = Fast
# LREAL = False
ISMEAR = 1
SIGMA = 0.2
# LASPH = True



# geometry optimization
EDIFFGG = -0.03
IBRION = -1
# ISIF = 2 # to check
NSW = 0
# ISPIN = 1
# LORBIT = 11
```

## Distribute Files to Folders

```
cp POSCAR POTCAR INCAR ./ktest
cp POSCAR POTCAR INCAR ./encut
cp runKPOINTS.py extractDataKPOINTS.py job ./ktest
cp runENCUT.py extractDataENCUT.py job ./encut
```

At this point you are ready to run the calculations! I hope this was helpful in getting you started.
