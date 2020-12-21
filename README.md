# SNfold
### a sequential niche multimodal conformation sampling algorithm for protein structure prediction



**Developer:**   
                Yuhao Xia  
                College of Information Engineering  
                Zhejiang University of Technology, Hangzhou 310023, China  
		Email: xiayh@zjut.edu.cn  

**Contact:**  
                Guijun Zhang, Prof  
                College of Information Engineering  
                Zhejiang University of Technology, Hangzhou 310023, China  
                Email: zgj@zjut.edu.cn  

## 1. INSTALLATION
Binaries for Linux 64 bit system has been included in the package. The Linux binary was compiled using GCC 5.4.0. Users need to have these versions of GCC compilers when using binaries.

Please Follow the below steps to install and configure SNfold:

- Download Rosetta3.10 source package from https://www.rosettacommons.org/software/.
(where `$ROSETTA3`=path-to-Rosetta)

- Copy and paste ``"ClassicAbinitio.cc"`` and ``"ClassicAbinitio.hh"`` from ``"src/"`` folder in SNfold package to ``"$ROSETTA3/main/source/src/protocols/abinitio/"`` folder in Rosetta.

- Copy and paste ``"MonteCarlo.cc"`` and ``"MonteCarlo.hh"`` from ``"src/"`` folder in SNfold package to ``"$ROSETTA3/main/source/src/protocols/moves/"`` folder in Rosetta.

- Compile SNfold source code using the following commands:

```
 $ cd $ROSETTA3/main/source/
 $ ./scons.py AbinitioRelax -j<NumOfJobs> mode=release bin
```

## 2. INPUT
SNfold requires four files to generate models:

	-fasta				: fasta file
	-alldist.txt			: distance map file
	-aat000_03_05.200_v1_3		: fragment library with fragment lenth 3
	-aat000_09_05.200_v1_3		: fragment library with fragment lenth 9

## 3. RUN
Please follow the below steps to run SNfold:

- Go to the ``"example/"`` folder of SNfold.

- Run SNfold with the following command:

```
   $ $ROSETTA3/main/source/bin/AbinitioRelax.default.linuxgccrelease @flags
```

- Five models are generated in the ``"output_files/"`` folder.

- ``"score.fsc"``, ``"S_00000001.pdb"`` and ``"default.out"`` can be deleted.

## 4. OUTPUT
Output files of SNfold are stored in the ``"example/output_files/"`` folder, including five predicted models (model_X.pdb).

	model_1.pdb
	model_2.pdb
	model_3.pdb
	model_4.pdb
	model_5.pdb

- ``"TMscore"`` in the ``"output_files/"`` folder can be used to calculate the accuracy of predicted models using the following commands:

```
 $ ./TMscore model_X.pdb ../input_files/native.pdb
```

## 5. DISCLAIMER
The executable software and the source code of SNfold is distributed free of charge 
as it is to any non-commercial users. The authors hold no liabilities to the performance 
of the program.
