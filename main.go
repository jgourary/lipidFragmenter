package main

import (
	"bufio"
	"errors"
	"fmt"
	"log"
	"os"
	filepath2 "path/filepath"
	"strconv"
	"strings"
	"time"
)

// Hydrogen bond length in angstroms. Used to add hydrogens to molecule fragments after cutting the molecule
const hydrogenBondDistance float64 = 1.10

func main() {
	t1 := time.Now()
	args := os.Args
	filePath := args[1]

	fmt.Println("\nLoading lipid: " + filePath)
	atoms := loadLipid(filePath)
	fmt.Println("Lipid loaded.")

	fmt.Println("\nAssigning heteroatoms & neighbors to groups...")
	createFunctionalGroups(atoms)
	fmt.Println("Assignment complete.")

	fmt.Println("\nAssigning alkane carbons to groups...")
	mergeAlkanes(atoms)
	fmt.Println("Alkane carbons assigned.")

	fmt.Println("\nAssigning hydrogens to groups...")
	mergeHydrogens(atoms)
	fmt.Println("Hydrogens assigned.")

	fmt.Println("\nFragmenting molecule...")
	separateGroups(atoms)
	fmt.Println("Molecule fragmented.")

	fmt.Println("\nWriting fragments to disk...")
	groups := getGroups(atoms)
	for i, group := range groups {
		groupMap, groupSlice := groupMap(atoms, group)
		writeFragment(atoms, groupMap, groupSlice, filePath, i)
	}
	fmt.Println("Fragments written.")
	dt := time.Since(t1).String()
	fmt.Println("\nProcess completed in " + dt)
}

func loadLipid(filePath string) map[int]*atom {

	// Create structure to store atoms
	atoms := make(map[int]*atom)

	// open file
	file, err := os.Open(filePath)
	if err != nil {
		fmt.Println("Failed to open key file: " + filePath)
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

func writeFragment(atoms map[int]*atom, groupMap map[int]int, groupSlice []int, filepath string, k int) {

	dir, filename := filepath2.Split(filepath)
	baseFileName := strings.Split(filename,".")[0]

	thisPath := filepath2.Join(dir,baseFileName + "_fragment_" + strconv.Itoa(k) + ".xyz")
	thisFile, err := os.Create(thisPath)
	if err != nil {
		fmt.Println("Failed to create new fragment file: " + thisPath)
		log.Fatal(err)
	}

	// write header
	_, err = thisFile.WriteString(strconv.Itoa(len(groupMap)) + "\t " + baseFileName + " Fragment " + strconv.Itoa(k) + "\n")
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

type atom struct {
	element string
	atomType int
	bondedAtoms []int

	pos []float64

	// used for union-find
	parent int
	treeSize int
	isInFuncGroup bool
}

func union(atoms map[int]*atom, atom1 int, atom2 int) {
	root1 := root(atoms, atom1)
	root2 := root(atoms, atom2)
	if root1 != root2 {
		if atoms[root1].treeSize < atoms[root2].treeSize {
			atoms[root1].parent = root2
			atoms[root2].treeSize += atoms[root1].treeSize
		} else {
			atoms[root2].parent = root1
			atoms[root1].treeSize += atoms[atom1].treeSize
		}
	}
}

func connected(atoms map[int]*atom, atom1 int, atom2 int) bool {
	root1 := root(atoms, atom1)
	root2 := root(atoms, atom2)
	if root1 != root2 {
		return false
	} else {
		return true
	}
}

func validate(atoms map[int]*atom, atom1 int) bool {
	if _, ok := atoms[atom1]; ok {
		return true
	} else {
		newErr := errors.New("Atom key " + strconv.Itoa(atom1) + " was not found in map")
		log.Fatal(newErr)
	}
	return false
}

func root(atoms map[int]*atom, atom1 int) int {
	validate(atoms, atom1)
	// store array of visited atoms for path compression afterwards
	var visitedAtoms []int

	// check if atom's parent is equal to atom's parent's parent (i.e. we have reached the top of the tree)
	for atoms[atom1].parent != atoms[atoms[atom1].parent].parent {
		// if not set atom's parent to atom's parent's parent
		atoms[atom1].parent = atoms[atoms[atom1].parent].parent
		// save atom's parent to list of visited atoms to path compress afterwards
		visitedAtoms = append(visitedAtoms,atoms[atom1].parent)
	}

	// compress path
	for _,visitedAtom := range visitedAtoms {
		atoms[visitedAtom].parent = atoms[atom1].parent
	}

	return atoms[atom1].parent
}
