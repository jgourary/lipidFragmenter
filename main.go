package main

import (
	"fmt"
	"io/ioutil"
	"log"
	"os"
	filepath2 "path/filepath"
)

// OPENBABEL exe locations
const obabel string = "C:\\Program Files\\OpenBabel-3.1.1\\obabel.exe"
const obminimize string = "C:\\Program Files\\OpenBabel-3.1.1\\obminimize.exe"


// Carbon-Carbon bond length in angstroms. Used to add hydrogens to molecule fragments after cutting the molecule
const carbonCarbonBondDistance float64 = 1.54
// Hydrogen-Carbon bond length in angstroms. Used to add hydrogens to molecule fragments after cutting the molecule
const hydrogenCarbonBondDistance float64 = 1.10

// How many go routines to launch at once
const goRoutinesBatchSize int = 128
// the number of monomers to be incorporated
const topFragsNum int = 100
// How many times a dimer must appear to be incorporated
const dimerCutoff int = 10

// Program Control Variables for main() function
const libraryGenMode bool = true

	const readDatabase bool = false
	const smi2txyz bool = false
	const fragment bool = false
	const xyz2smi bool = false
	const procFreqs bool = false
	const generateLibraries bool = true
	const singleFragLimit int = 100
	const doubleFragLimit int = 25

const bilayerConversionMode bool = false

	const generateCodeDict bool = true
	const reviseBilayers bool = false

// Program begins here
func main() {
	filePath := "C:\\Users\\jtgou\\lipids2\\structures.sdf"

	// filePath := "/home/jtg2769/lipids/LMSD_20191002.sdf"
	// dir := "/home/jtg2769/lipids"

	// The primary program output directory where all the below subdirectories are held
	dir := "C:\\Users\\jtgou\\lipids2"
	// The below are subdirectories of the above for specific purposes

	// Stores the raw lipid molecules from the database in TXYZ form
	moleculesDir := filepath2.Join(dir,"molecules")
	// Stores the single and double fragments generated from the fragmentation of the database molecules
	singleFragmentsDir := filepath2.Join(dir,"single_fragments")
	doubleFragmentsDir := filepath2.Join(dir,"double_fragments")
	dimersDir := filepath2.Join(dir,"dimers")
	uniqueSingleFragsDir := filepath2.Join(dir,"unique_single_fragments")
	uniqueDoubleFragsDir := filepath2.Join(dir,"unique_double_fragments")

	library := filepath2.Join(dir,"fragment_library")
	librarySFcatalog := filepath2.Join(dir,"fragment_library","single_fragments.txt")
	libraryDFcatalog := filepath2.Join(dir,"fragment_library","double_fragments.txt")
	librarySFDir := filepath2.Join(library,"single_fragments")
	libraryDFDir := filepath2.Join(library,"double_fragments")

	// uniqueDimersDir := filepath2.Join(dir,"unique_dimers")
	os.Mkdir(dir,0755)


	if libraryGenMode {
		if readDatabase {
			fmt.Println("Reading database and converting contents to SMI molecule files...")
			readLipidMaps(filePath, moleculesDir)
		}
		if smi2txyz {
			fmt.Println("Converting SMI molecule files to TXYZ molecule files...")
			obabelConversion2(moleculesDir, ".smi", ".txyz", "add", false)
		}
		if fragment {
			fmt.Println("Dividing TXYZ molecule files into TXYZ fragments...")
			fragmentManager(moleculesDir, singleFragmentsDir, doubleFragmentsDir, dimersDir)
		}
		if xyz2smi {
			fmt.Println("Converting TXYZ fragments to CAN (unique SMILES) fragments...")
			obabelConversion2(singleFragmentsDir, ".txyz", ".sdf", "remove", false)
			obabelConversion2(singleFragmentsDir, ".sdf", ".can", "no", false)
			// obabelConversion2(doubleFragmentsDir, ".txyz", ".can", "remove", false)
		}
		if procFreqs {
			fmt.Println("Counting single and double fragment occurrences from CAN fragments...")
			fragSelector(dir, singleFragmentsDir, doubleFragmentsDir, uniqueSingleFragsDir, uniqueDoubleFragsDir)
		}
		if generateLibraries {
			fmt.Println("Generating library of most common single fragments TXYZs")
			inPath := filepath2.Join(dir, "top_single_fragments.txt")
			generateLibrary(inPath, librarySFcatalog, librarySFDir, singleFragLimit)

			// make SDFs for POLTYPE
			obabelConversion(librarySFDir, ".txyz", ".sdf", "add", true, false)
			createPoltypeINIs(librarySFDir)


			fmt.Println("Generating library of most common double fragments TXYZs")
			inPath = filepath2.Join(dir, "top_double_fragments.txt")
			generateLibrary(inPath, libraryDFcatalog, libraryDFDir, doubleFragLimit)

			// make SDFs for POLTYPE
			obabelConversion(libraryDFDir, ".txyz", ".sdf", "add", true, false)
			createPoltypeINIs(libraryDFDir)
		}
	} else if bilayerConversionMode {
		dir := "C:\\Users\\jtgou\\OneDrive\\Documents\\UT_Austin\\ren_lab\\lipids\\pren\\molecules"
		outDir := "C:\\Users\\jtgou\\OneDrive\\Documents\\UT_Austin\\ren_lab\\lipids\\pren\\atomCodeDict"
		outName := "atomCodeDict.txt"
		if generateCodeDict {
			generateAtomCodeDictFile(dir, outDir, outName)
		}

		if reviseBilayers {
			bilayerInDir := "C:\\Users\\jtgou\\OneDrive\\Documents\\UT_Austin\\ren_lab\\lipids\\pren\\molecules"
			bilayerOutDir := "C:\\Users\\jtgou\\OneDrive\\Documents\\UT_Austin\\ren_lab\\lipids\\pren\\bilayers"
			atomCodeFile := filepath2.Join(outDir, outName)
			processBilayers(bilayerInDir, atomCodeFile, bilayerOutDir)
		}
	}

}

func createPoltypeINIs(dir string) {
	fileInfo, err := ioutil.ReadDir(dir)
	if err != nil {
		fmt.Println("failed to read directory: " + dir)
		log.Fatal(err)
	}
	for i := 0; i < len(fileInfo); i++ {

		if filepath2.Ext(fileInfo[i].Name()) == ".sdf" {
			sdfName := fileInfo[i].Name()
			iniName := "poltype.ini"
			iniFile := filepath2.Join(dir, iniName)
			//fmt.Println(iniFile)
			createPoltypeINI(iniFile, sdfName)
		} else if fileInfo[i].IsDir() {
			createPoltypeINIs(filepath2.Join(dir, fileInfo[i].Name()))
		}
	}
}

func createPoltypeINI(thisPath string, sdfName string) {
	thisFile, err := os.Create(thisPath)
	if err != nil {
		fmt.Println("Failed to create new fragment file: " + thisPath)
		log.Fatal(err)
	}

	// write header
	_, err = thisFile.WriteString("structure=" + sdfName + "\n")
	_, err = thisFile.WriteString("numproc=4\n")
	_, err = thisFile.WriteString("maxmem=20GB\n")
	_, err = thisFile.WriteString("maxdisk=100GB\n")
	_, err = thisFile.WriteString("externalapi=RenLabCluster\n")
	_, err = thisFile.WriteString("username=jtg2769\n")
}
