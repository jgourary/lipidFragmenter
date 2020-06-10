package main

import (
	"fmt"
	"os"
	filepath2 "path/filepath"
)

// Hydrogen bond length in angstroms. Used to add hydrogens to molecule fragments after cutting the molecule
const carbonCarbonBondDistance float64 = 1.54
const hydrogenCarbonBondDistance float64 = 1.10

const goRoutinesBatchSize int = 128

const readDatabase bool = false
const sdf2xyz bool = false
const fragment bool = false
const xyz2smi bool = false
const procSMIs bool = true

func main() {
	filePath := "C:\\Users\\jtgou\\OneDrive\\Documents\\UT_Austin\\ren_lab\\lipids\\lipidmaps\\LMSD_20191002.sdf"
	dir := "C:\\Users\\jtgou\\OneDrive\\Documents\\UT_Austin\\ren_lab\\lipids\\lipidmaps"
	moleculesDir := filepath2.Join(dir,"molecules")
	fragmentsDir := filepath2.Join(dir,"fragments")
	os.Mkdir(dir,0755)

	fmt.Println("Reading database and converting contents to SDF files...")
	if readDatabase { readLipidMaps(filePath, moleculesDir) }
	fmt.Println("Converting SDF files to TXYZ files...")
	if sdf2xyz { obabelConversion(moleculesDir, ".sdf", ".txyz", true) }
	fmt.Println("Dividing TXYZ files into TXYZ fragments...")
	if fragment { fragmentManager(moleculesDir, fragmentsDir) }
	fmt.Println("Converting XYZ fragments to SMI fragments...")
	if xyz2smi { obabelConversion(fragmentsDir, ".txyz", ".can", false) }
	fmt.Println("Counting fragment occurrences...")
	if procSMIs { smiProcessor(dir, fragmentsDir) }


}

