# CGLFold
### a contact-assisted de novo protein structure prediction using global exploration and loop perturbation sampling algorithm

**Developer:**   
                Jun Liu  
                College of Information Engineering  
                University of Zhejiang University of Technology, Zhejiang  
                Email: junl@zjut.edu.cn  
		
**Contact:**  
                Guijun Zhang, Prof  
                College of Information Engineering  
                University of Zhejiang University of Technology, Zhejiang  
                Email: zgj@zjut.edu.cn  

## 1. INSTALLATION
Binaries for Linux 64 bit system has been included in the package. The Linux binary was compiled using GCC 5.4.0. Users need to have these versions of GCC compilers when using binaries.

Please Follow the below steps to install and configure CGLFold:

- Download Rosetta3.10 source package from https://www.rosettacommons.org/software/ 
and extract it to ``"~/"`` directory.

- Copy and paste source code of ``"ClassicAbinitio.cc"``, ``"ClassicAbinitio.hh"``,   ``"LJAngleRotation.cc"``, ``"LJAngleRotation.hh"``, and ``"LJAngleRotation.fwd.hh"`` from ``"src/"`` folder in CGLFold package to ``"~/Rosetta/main/source/src/protocols/abinitio/"`` folder in Rosetta. Copy and paste configuration file ``"protocols_b_6.src.settings"`` from ``"src/"`` folder in CGLFold package to ``"~/Rosetta/main/source/src/"`` folder in Rosetta.

- Compile Rosetta source code using the following commands:

```
 $ cd ~/Rosetta/main/source/
 $ ./scons.py -j<NumOfJobs> mode=release bin
```

- If you want to recompile CGLFold source code, use the following commands:

```
 $ cd ~/CONFold/
 $ g++ -o bin/CGLFold src/CGLFold.cpp
```
## 2. INPUT
CGLFold requires four files to generate models:

	-f	fasta			: fasta file
	-c	cmap			: contact map file
	-frag3	3mer_fragment_library	: fragment library with fragment lenth 3
	-frag9	9mer_fragment_library	: fragment library with fragment lenth 9

## 3. OUTPUT
Output files of CGLFold are stored in the ``"example/output_files/"`` folder, including five predicted model (model_X.pdb) and the filtered contact map.

	model_1.pdb		first cluster centroid structure obtained by SPICKER
	model_2.pdb		the structure closest to first cluster centroid
	model_3.pdb		the structure with lowest Cscore in the last generation
	model_4.pdb		the structure with lowest energy in the last generation
	model_5.pdb		secondary cluster centroid structure obtained by SPICKER
	filtered_camp.txt	used to build the contact-based selection model

## 4. EXAMPLE
Please follow the below steps to run CGLFold:

- Go to the ``"example/"`` folder of CGLFold.
  
- Run CGLFold with the following command:
  
```
   $ ./../bin/CGLFold -f ./input_files/fasta.txt -c input_files/contact.txt -frag3 input_files/3mer_fragment_library -frag9 input_files/9mer_fragment_library
```

- Five models and the filtered contact map are generated in the ``"output_files/"`` folder.

## 5. DISCLAIMER
The executable software and the source code of CGLFold is distributed free of charge 
as it is to any non-commercial users. The authors hold no liabilities to the performance 
of the program.
