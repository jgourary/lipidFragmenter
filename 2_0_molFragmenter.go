package main

import (
	"bufio"
	"errors"
	"fmt"
	"io/ioutil"
	"log"
	"os"
	filepath2 "path/filepath"
	"strconv"
	"strings"
	"sync"
)

func fragmentManager(moleculesDir string, singleFragmentsDir string, doubleFragmentsDir string, dimersDir string) {

	// remove existing fragments
	singleFragsExist, _ := exists(singleFragmentsDir)
	if singleFragsExist {
		_ = os.RemoveAll(singleFragmentsDir)
	}
	doubleFragsExist, _ := exists(doubleFragmentsDir)
	if doubleFragsExist {
		_ = os.RemoveAll(doubleFragmentsDir)
	}
	dimersExist, _ := exists(dimersDir)
	if dimersExist {
		_ = os.RemoveAll(dimersDir)
	}

	fileInfo, err := ioutil.ReadDir(moleculesDir)
	if err != nil {
		fmt.Println("failed to read directory: " + moleculesDir)
		log.Fatal(err)
	}

	maxFrame := min(goRoutinesBatchSize,len(fileInfo)-1)
	frame := []int{0,maxFrame}

	for frame[0] < len(fileInfo) {

		wg := sync.WaitGroup{}

		for i := frame[0]; i <= frame[1]; i++ {
			if filepath2.Ext(fileInfo[i].Name()) == ".txyz" {
				molPath := filepath2.Join(moleculesDir, fileInfo[i].Name())

				wg.Add(1)
				go fragmentMoleculeShellFunc(molPath, singleFragmentsDir, doubleFragmentsDir, dimersDir, &wg)
				// fragmentMoleculeShellFunc(molPath, singleFragmentsDir, doubleFragmentsDir, dimersDir, &wg)

			}
		}


		frame[0] += goRoutinesBatchSize
		frame[1] += goRoutinesBatchSize
		frame[1] = min(frame[1], len(fileInfo)-1)
		wg.Wait()
	}

}

func fragmentMoleculeShellFunc(filePath string, singleFragmentsDir string, doubleFragmentsDir string, dimersDir string, wg *sync.WaitGroup) {

	lipidName, atoms := loadLipid(filePath)
	fragmentMolecule(atoms, lipidName, singleFragmentsDir, doubleFragmentsDir, dimersDir)

	wg.Done()

}

func fragmentMolecule(atoms map[int]*atom, lipidName string, singleFragmentsDir string, doubleFragmentsDir string, dimersDir string) {

	// fmt.Println("Fragmenting lipid " + lipidName)
	//fmt.Println("\nLooking for bridge bonds...")
	bridge(atoms)
	// fmt.Println("Bridge bonds and acyclic atoms identified.")

	//fmt.Println("\nAssigning heteroatoms & neighbors to groups...")
	createFunctionalGroups(atoms)
	// fmt.Println("Assignment complete.")

	//fmt.Println("\nAssigning alkane carbons to groups...")
	mergeAlkanes(atoms)
	// fmt.Println("Alkane carbons assigned.")

	//fmt.Println("\nAssigning hydrogens to groups...")
	mergeHydrogens(atoms)
	// mt.Println("Hydrogens assigned.")

	//fmt.Println("\nIdentifying borders between fragments...")
	borderBonds := getFragmentBorderBonds(atoms)
	// fmt.Println("Border Identified.")

	//fmt.Println("\nFinding single fragments...")
	singleFrags, _ := getSingleFragments(atoms,borderBonds)
	// fmt.Println("Single Fragments Found.")

	//fmt.Println("\nFinding double fragments...")
	doubleFrags, dimers, _ := getDoubleFragments(atoms,borderBonds)
	// doubleFrags, doubleFragsIsHydrocarbon := getDoubleFragments(atoms,borderBonds)
	// fmt.Println("BDouble Fragments found.")

	//fmt.Println("\nWriting single fragments to disk...")
	// fmt.Println(len(singleFrags))
	for i := 0; i < len(singleFrags); i++ {
		fragRoot := root(singleFrags[i], 1)
		fragName := lipidName + "_single_" + strconv.Itoa(fragRoot)
		writeFragment(singleFrags[i], singleFragmentsDir, fragName)
	}

	//fmt.Println("\nWriting double fragments to disk...")
	for i := 0; i < len(doubleFrags); i++ {
		fragRoot1 := root(atoms, borderBonds[i][0])
		fragRoot2 := root(atoms, borderBonds[i][1])
		fragName := lipidName + "_double_" + strconv.Itoa(fragRoot1) + "_" + strconv.Itoa(fragRoot2)
		writeFragment(doubleFrags[i], doubleFragmentsDir, fragName)
	}

	//fmt.Println("\nWriting double fragments to disk...")
	for i := 0; i < len(dimers); i++ {
		fragRoot1 := root(atoms, borderBonds[i][0])
		fragRoot2 := root(atoms, borderBonds[i][1])
		fragName := lipidName + "_dimer_" + strconv.Itoa(fragRoot1) + "_" + strconv.Itoa(fragRoot2)
		writeFragment(dimers[i], dimersDir, fragName)
	}
}

func loadLipid(filePath string) (string, map[int]*atom) {

	// Create structure to store atoms
	atoms := make(map[int]*atom)

	// open file
	file, err := os.Open(filePath)
	if err != nil {
		fmt.Println("Failed to open molecule file: " + filePath)
		log.Fatal(err)
	}
	lipidName := strings.Split(filepath2.Base(filePath),".")[0]


	// Initialize scanner
	scanner := bufio.NewScanner(file)
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

			// create new atom
			var newAtom atom

			// get number of atom from file
			atomNum, err := strconv.Atoi(tokens[0])
			if err != nil {
				newErr := errors.New("Failed to convert token in position 0 on line " + strconv.Itoa(i) + " to an integer")
				log.Fatal(newErr)
			}

			// assign parent and size for union-find algorithm
			newAtom.parent = atomNum
			newAtom.treeSize = 1

			// assign element
			newAtom.element = tokens[1]

			// assign positions
			pos := make([]float64,3)
			for j := 2; j < 5; j++ {
				pos[j-2], err = strconv.ParseFloat(tokens[j],64)
				if err != nil {
					newErr := errors.New("Failed to convert token in position 0 on line " + strconv.Itoa(j) + " to a float64")
					log.Fatal(newErr)
				}
			}
			newAtom.pos = pos

			// assign atomType from file
			newAtom.atomType, err = strconv.Atoi(tokens[5])
			if err != nil {
				newErr := errors.New("Failed to convert token in position 5 on line " + strconv.Itoa(i) + " to an integer")
				log.Fatal(newErr)
			}

			// assign bonds from file
			bonds := make([]int,len(tokens)-6)
			for j := 6; j < len(tokens); j++ {
				bonds[j-6], err = strconv.Atoi(tokens[j])
				if err != nil {
					newErr := errors.New("Failed to convert token in position " + strconv.Itoa(j) + " on line " + strconv.Itoa(i) + " to an integer")
					log.Fatal(newErr)
				}
			}
			newAtom.bondedAtoms = bonds

			// assign parameters for bridge finding alg
			newAtom.discTime = 1e10
			newAtom.minDiscTime = 1e10
			newAtom.parentBF = -1
			newAtom.visited = false
			newAtom.isCyclic = false

			// add atom to map
			atoms[atomNum] = &newAtom

		} else {
			fmt.Println("Warning: line " + strconv.Itoa(i) + " has insufficient tokens. Program is skipping this " +
				"line when reading your input file.")
		}
		i++
	}

	return lipidName, atoms
}

func getFragmentCharge(atoms map[int]*atom) int {
	charge := 0
	for _, atom := range atoms {
		if atom.element == "N" {
			if len(atom.bondedAtoms) == 4 {
				charge++
			}
		} else if atom.element == "P" {
			incrementCharge := false
			for _, atomID := range atom.bondedAtoms {
				if atoms[atomID].element == "O" &&  len(atoms[atomID].bondedAtoms) == 1 {
					incrementCharge = true
				}
			}
			if incrementCharge {
				charge--
			}
		}
	}
	return charge
}

func writeFragment(atoms map[int]*atom, fragSubDir string, fragName string) {
	os.MkdirAll(fragSubDir, 0755)
	thisPath := filepath2.Join(fragSubDir, fragName + ".txyz")
	// fmt.Println(thisPath)
	thisFile, err := os.Create(thisPath)
	if err != nil {
		fmt.Println("Failed to create new fragment file: " + thisPath)
		log.Fatal(err)
	}

	// write header
	_, err = thisFile.WriteString(strconv.Itoa(len(atoms)) + "\t Fragment " + fragName + "charge=" + strconv.Itoa(getFragmentCharge(atoms)) + "\n")
	if err != nil {
		fmt.Println("Failed to write header line to key: " + thisPath)
		log.Fatal(err)
	}

	// write body
	for i := 1; i <= len(atoms); i++ {
		line := strconv.Itoa(i) + "\t" + atoms[i].element + "\t" + fmt.Sprintf("%.6f",atoms[i].pos[0]) + "\t" +
			fmt.Sprintf("%.6f",atoms[i].pos[1]) + "\t" + fmt.Sprintf("%.6f",atoms[i].pos[2]) + "\t" +
			strconv.Itoa(atoms[i].atomType)
		for _, bondedAtom := range atoms[i].bondedAtoms {
			line += "\t" + strconv.Itoa(bondedAtom)
		}

		_, err = thisFile.WriteString(line + "\n")
		if err != nil {
			fmt.Println("Failed to write header line to key: " + thisPath)
			log.Fatal(err)
		}
	}

}

