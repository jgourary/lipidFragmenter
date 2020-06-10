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

func fragmentManager(moleculesDir string, fragmentsDir string) {

	fileInfo, err := ioutil.ReadDir(moleculesDir)
	if err != nil {
		fmt.Println("failed to read directory: " + moleculesDir)
		log.Fatal(err)
	}

	frame := []int{0,goRoutinesBatchSize-1}
	for frame[0] < len(fileInfo) {

		wg := sync.WaitGroup{}

		for i := frame[0]; i < frame[1]; i++ {
			// if item is a Dir (as it should be unless the end user tampered with the directory manually...)
			if fileInfo[i].IsDir() {
				molSubDir := filepath2.Join(moleculesDir,fileInfo[i].Name())
				fragSubDir := filepath2.Join(fragmentsDir,fileInfo[i].Name())

				molSubDirInfo, err := ioutil.ReadDir(molSubDir)
				if err != nil {
					fmt.Println("failed to read directory: " + molSubDir)
					log.Fatal(err)
				}
				for j := 0; j < len(molSubDirInfo); j++ {
					fileName := molSubDirInfo[j].Name()
					if filepath2.Ext(fileName) == ".txyz" {
						molPath := filepath2.Join(molSubDir, fileName)
						go fragmentMolecule(molPath, fragSubDir, &wg)
						// Add one to wait group
						wg.Add(1)
					}
				}
			}
		}

		wg.Wait()

		frame[0] += goRoutinesBatchSize
		frame[1] += goRoutinesBatchSize
		frame[1] = min(frame[1], len(fileInfo))
	}

}

func fragmentMolecule(filePath string, fragSubDir string, wg *sync.WaitGroup) {

	// remove existing fragments
	fragsExist, _ := exists(fragSubDir)
	if fragsExist {
		os.RemoveAll(fragSubDir)
	}

	//fmt.Println("\nLoading lipid: " + filePath)
	atoms := loadLipid(filePath)
	//fmt.Println("Lipid loaded.")

	//fmt.Println("\nLooking for bridge bonds...")
	bridge(atoms)
	//fmt.Println("Bridge bonds and acyclic atoms identified.")

	//fmt.Println("\nAssigning heteroatoms & neighbors to groups...")
	createFunctionalGroups(atoms)
	//fmt.Println("Assignment complete.")

	//fmt.Println("\nAssigning alkane carbons to groups...")
	mergeAlkanes(atoms)
	//fmt.Println("Alkane carbons assigned.")

	//fmt.Println("\nAssigning hydrogens to groups...")
	mergeHydrogens(atoms)
	//fmt.Println("Hydrogens assigned.")

	//fmt.Println("\nFragmenting molecule...")
	separateGroups(atoms)
	//fmt.Println("Molecule fragmented.")

	//fmt.Println("\nWriting fragments to disk...")
	groups := getGroups(atoms)
	j := 0
	for _, group := range groups {
		groupMap, groupSlice := groupMap(atoms, group)
		writeFragment(atoms, groupMap, groupSlice, fragSubDir, j)
		j++
	}
	//fmt.Println("Fragments written.")

	wg.Done()

}

func loadLipid(filePath string) map[int]*atom {

	// Create structure to store atoms
	atoms := make(map[int]*atom)

	// open file
	file, err := os.Open(filePath)
	if err != nil {
		fmt.Println("Failed to open molecule file: " + filePath)
		log.Fatal(err)
	}
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

	return atoms
}

func writeFragment(atoms map[int]*atom, groupMap map[int]int, groupSlice []int, fragSubDir string, k int) {
	os.MkdirAll(fragSubDir, 0755)
	thisPath := filepath2.Join(fragSubDir, "fragment_" + strconv.Itoa(k) + ".txyz")
	// fmt.Println(thisPath)
	thisFile, err := os.Create(thisPath)
	if err != nil {
		fmt.Println("Failed to create new fragment file: " + thisPath)
		log.Fatal(err)
	}

	// write header
	_, err = thisFile.WriteString(strconv.Itoa(len(groupMap)) + "\t Fragment " + strconv.Itoa(k) + "\n")
	if err != nil {
		fmt.Println("Failed to write header line to key: " + thisPath)
		log.Fatal(err)
	}

	// write body
	for _, atomID := range groupSlice {
		line := strconv.Itoa(groupMap[atomID]) + "\t" + atoms[atomID].element + "\t" + fmt.Sprintf("%.6f",atoms[atomID].pos[0]) + "\t" +
			fmt.Sprintf("%.6f",atoms[atomID].pos[1]) + "\t" + fmt.Sprintf("%.6f",atoms[atomID].pos[2]) + "\t" +
			strconv.Itoa(atoms[atomID].atomType)
		for _, bondedAtom := range atoms[atomID].bondedAtoms {
			line += "\t" + strconv.Itoa(groupMap[bondedAtom])
		}

		_, err = thisFile.WriteString(line + "\n")
		if err != nil {
			fmt.Println("Failed to write header line to key: " + thisPath)
			log.Fatal(err)
		}
	}

}

