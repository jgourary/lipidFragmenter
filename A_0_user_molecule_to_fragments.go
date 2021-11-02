package main

import (
	"bufio"
	"fmt"
	"io/ioutil"
	"log"
	"os"
	filepath2 "path/filepath"
	"strings"
	"sync"
)

func conversionManager(inFilePath string, outDir string) {

	lipidName, atoms := loadLipid(inFilePath)

	singleFragmentsDir := filepath2.Join(outDir,"single_fragments")
	doubleFragmentsDir := filepath2.Join(outDir,"double_fragments")
	dimerFragmentsDir := filepath2.Join(outDir,"dimers")

	fragmentMolecule(atoms, lipidName, singleFragmentsDir, doubleFragmentsDir, dimerFragmentsDir)

	obabelConversion(singleFragmentsDir, ".txyz", ".can", "no", false, true)
	obabelConversion(doubleFragmentsDir, ".txyz", ".can", "no", false, true)
	// obabelConversion(dimerFragmentsDir, ".txyz", ".can", false, false)

	fragDatabaseDir := "placeholder"

	singleFragDatabase, _ := loadFragmentDatabase(fragDatabaseDir)
	// singleFragDatabase, doubleFragDatabase := loadFragmentDatabase(fragDatabaseDir)

	// Read in all files in dir
	fileInfo, err := ioutil.ReadDir(singleFragmentsDir)
	if err != nil {
		fmt.Println("failed to read directory: " + singleFragmentsDir)
		log.Fatal(err)
	}

	// create out file
	thisPath := filepath2.Join(outDir, lipidName + ".out")
	// fmt.Println(thisPath)
	thisFile, err := os.Create(thisPath)
	if err != nil {
		fmt.Println("Failed to create new fragment file: " + thisPath)
		log.Fatal(err)
	}
	_, _ = thisFile.WriteString("Lipid Fragmenter Output - " + lipidName + "\n")


	for i := 0; i < len(fileInfo); i++ {
		if filepath2.Ext(fileInfo[i].Name()) == ".can" {

			canFilePath := filepath2.Join(singleFragmentsDir, fileInfo[i].Name())

			smiString := getSMIString(canFilePath)

			// record which fragment was matched in log file
			if path, ok := singleFragDatabase[smiString]; ok {
				_, _ = thisFile.WriteString(smiString + " " + path)
			} else {
				_, _ = thisFile.WriteString(smiString + " " + "No matching fragment found")
			}
		}
	}
}

func getAtomIDsToAtomTypesMap(atoms map[int]*atom, atomCodeDict map[string]int) map[int]int {
	atomIDtoType := make(map[int]int)
	wg := sync.WaitGroup{}
	for atomID := range atoms {
		wg.Add(1)
		thisID := atomID
		go func(wg *sync.WaitGroup) {
			atomCode := getAtomCode(atoms, thisID)
			atomType := atomCodeDict[atomCode]
			atomIDtoType[thisID] = atomType
			wg.Done()
		}(&wg)
	}
	wg.Wait()
	return atomIDtoType

}

func getAtomCode(atoms map[int]*atom, atomID int) string {

	atomPathways := make([]string, len(atoms[atomID].bondedAtoms))

	for i, bondedAtomID := range atoms[atomID].bondedAtoms {
		bondedAtomElement := atoms[bondedAtomID].element

		secondTierBondedElements := make([]string, len(atoms[bondedAtomID].bondedAtoms))
		for j, secondTierBondedAtomID := range atoms[bondedAtomID].bondedAtoms {
			secondTierBondedElements[j] = atoms[secondTierBondedAtomID].element
		}
		secondTierBondedElements = qsort4(secondTierBondedElements)

		atomPathways[i] = bondedAtomElement + "(" + strings.Join(secondTierBondedElements, "") + ")"

	}
	atomPathways = qsort4(atomPathways)

	atomIdentifier := atoms[atomID].element + "[" + strings.Join(atomPathways, "") + "]"

	return atomIdentifier
}

func loadFragmentDatabase(dir string) (map[string]string, map[string]string) {

	singleFragmentsDir := filepath2.Join(dir,"single_fragments")
	singleFragsMap := processDatabaseFolder(singleFragmentsDir)
	doubleFragmentsDir := filepath2.Join(dir,"double_fragments")
	doubleFragsMap := processDatabaseFolder(doubleFragmentsDir)

	return singleFragsMap, doubleFragsMap
}

func processDatabaseFolder(dir string) map[string]string {
	fragDirFileInfo, err := ioutil.ReadDir(dir)
	if err != nil {
		fmt.Println("failed to read directory: " + dir)
		log.Fatal(err)
	}

	smiles2filePath := make(map[string]string)



	fmt.Println("Identifying all non-alkane fragments...")

	for i := 0; i < len(fragDirFileInfo); i++ {
		fileName := fragDirFileInfo[i].Name()
		if filepath2.Ext(fileName) == ".can" {

			smiFilePath := filepath2.Join(dir, fileName)

			file, err := os.Open(smiFilePath)
			if err != nil {
				fmt.Println("Failed to open molecule file: " + smiFilePath)
				log.Fatal(err)
			}
			// Initialize scanner
			scanner := bufio.NewScanner(file)
			// ignore first line
			scanner.Scan()
			line := scanner.Text()
			fields := strings.Fields(line)
			if len(fields) > 1 {
				smiString := fields[0]
				filePath := fields[1]
				smiles2filePath[smiString] = filePath
			}

		}
	}

	return smiles2filePath
}