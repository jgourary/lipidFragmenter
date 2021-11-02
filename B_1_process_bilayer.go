package main

import (
	"bufio"
	"fmt"
	"io/ioutil"
	"log"
	"os"
	filepath2 "path/filepath"
	"strconv"
	"strings"
	"sync"
)

// Will convert assign correct atom types to all bilayers in a dir
func processBilayers(bilayerDir string, atomCodeFile string, outDir string) {
	// Read in all files in dir
	fileInfo, err := ioutil.ReadDir(bilayerDir)
	if err != nil {
		fmt.Println("failed to read directory: " + bilayerDir)
		log.Fatal(err)
	}
	wg := sync.WaitGroup{}
	for i := 0; i < len(fileInfo); i++ {
		if filepath2.Ext(fileInfo[i].Name()) == ".txyz" {
			wg.Add(1)
			go func(wg *sync.WaitGroup) {
				bilayerFile := filepath2.Join(bilayerDir, fileInfo[i].Name())
				processBilayer(bilayerFile, atomCodeFile, outDir)
				wg.Done()
			}(&wg)
		}
	}
	wg.Wait()
}

// Will convert assign correct atom types to one bilayer
func processBilayer(bilayerFile string, atomCodeFile string, outDir string) {

	fmt.Println("Loading next bilayer into memory...")
	// Load bilayer into memory
	bilayerName, bilayer := loadLipid(bilayerFile)
	fmt.Println("Finished loading bilayer " + bilayerName + " into memory.")
	fmt.Println()

	fmt.Println("Loading Atom Type Assignment Database...")
	// Load atom code to atom type database
	atomCodeDict := getAtomCodeDictFromFile(atomCodeFile)
	fmt.Println("Finished loading Atom Type Assignment Database.")
	fmt.Println()

	fmt.Println("Assigning atom types...")
	// Assign correct biotypes to all molecules using atom code dict
	atomIDsToTypesMap := getAtomIDsToAtomTypesMap(bilayer, atomCodeDict)
	fmt.Println("Finished assigning atom types.")
	fmt.Println()

	fmt.Println("Writing output to: " + outDir + " ...")
	rewriteBilayer(atomIDsToTypesMap, bilayerFile, outDir, bilayerName)
	fmt.Println("Finished writing output.")
	fmt.Println()

}

// Loads dictionary of atom codes from disk
func getAtomCodeDictFromFile(file string) map[string]int {
	// Create structure to store atoms
	atomCodeDict := make(map[string]int)

	// open file
	thisFile, err := os.Open(file)
	if err != nil {
		fmt.Println("Failed to open molecule file: " + file)
		log.Fatal(err)
	}

	// Initialize scanner
	scanner := bufio.NewScanner(thisFile)

	// iterate over all other lines
	for scanner.Scan() {
		// get next line
		line := scanner.Text()
		// split by whitespace
		tokens := strings.Fields(line)

		if len(tokens) > 2 {
			k := tokens[0]
			v := tokens[1]
			atomCodeDict[k], err = strconv.Atoi(v)
			if err != nil {
				fmt.Println("Warning: could not parse token " + tokens[1] + " as integer in file: " + file + " in line \"" + line + "\" at token position 2")
			}
		}
	}

	return atomCodeDict
}

// Rewrites file contents of bilayer with correct atom types
func rewriteBilayer(atomIDsToTypesMap map[int]int, bilayerFile string, outDir string, bilayerName string) {
	_ = os.MkdirAll(outDir, 0755)
	outPath := filepath2.Join(outDir, bilayerName + "_amoeba.txyz")
	// fmt.Println(thisPath)
	outFile, err := os.Create(outPath)
	if err != nil {
		fmt.Println("Failed to create new fragment file: " + outPath)
		log.Fatal(err)
	}
	// write header
	_, err = outFile.WriteString(strconv.Itoa(len(atomIDsToTypesMap)) + "\t System: " + bilayerName + "\n")
	if err != nil {
		fmt.Println("Failed to write header line to key: " + outPath)
		log.Fatal(err)
	}

	// open file
	inFile, err := os.Open(bilayerFile)
	if err != nil {
		fmt.Println("Failed to open molecule file: " + bilayerFile)
		log.Fatal(err)
	}

	// Initialize scanner
	scanner := bufio.NewScanner(inFile)
	// ignore first line
	scanner.Scan()
	// create line counter
	i := 1
	// iterate over all other lines
	for scanner.Scan() {
		// get next line
		line := scanner.Text()
		// split by whitespace
		tokens := strings.Fields(line)
		// check line length before proceeding
		if len(tokens) >= 6 {
			// get atomID
			atomID, err := strconv.Atoi(tokens[0])
			if err != nil {
				fmt.Println("Failed to read atom number on line " + strconv.Itoa(i) + "in file: " + bilayerFile)
				log.Fatal(err)
			}
			// get Atom Type from atomID
			newAtomType := atomIDsToTypesMap[atomID]
			// reassign Atom Type in line
			tokens[5] = strconv.Itoa(newAtomType)
		}
		// rewrite line
		replacementLine := ""
		for _, token := range tokens {
			replacementLine += token + "\t"
		}
		// write line to out file
		_, _ = outFile.WriteString(replacementLine + "\n")
		i++
	}

}

/*func getSeparateMolecules(atoms map[int]*atom) []map[int]*atom {
	roots := getRoots(atoms)
	moleculeSlice := make([]map[int]*atom, len(roots))

	// Parallelized iterate through all roots and create copy of molecule containing only atoms with that root
	wg := sync.WaitGroup{}
	for i := 0; i < len(roots); i++ {
		wg.Add(1)
		go func(wg *sync.WaitGroup) {
			// create a copyFile of the molecule
			moleculeCopy := copyMolecule(atoms)
			// delete all atoms not of that root
			deleteAtomsNotOfGroup(moleculeCopy, []int {roots[i]})
			// copy resulting molecule to the slice
			moleculeSlice[i] = moleculeCopy

			wg.Done()
		}(&wg)
	}
	wg.Wait()
	return moleculeSlice
}*/


