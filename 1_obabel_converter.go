package main

import (
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"os/exec"
	filepath2 "path/filepath"
	"strings"
	"sync"
)

// for file structure directory > subdirectory > file to be converted
func obabelConversion(directory string, ext1 string, ext2 string, addHydrogens string, addCoords bool, deleteOriginal bool) {
	// Read in all files in dir
	fileInfo, err := ioutil.ReadDir(directory)
	if err != nil {
		fmt.Println("failed to read directory: " + directory)
		log.Fatal(err)
	}


	// Iterate through all items in directory
	maxFrame := min(goRoutinesBatchSize,len(fileInfo)-1)
	frame := []int{0,maxFrame}

	for frame[0] < len(fileInfo) {

		wg := sync.WaitGroup{}

		for i := frame[0]; i <= frame[1]; i++ {
			// if item is a Dir - go through subdirs
			if fileInfo[i].IsDir() {
				subFileInfo, err := ioutil.ReadDir(filepath2.Join(directory, fileInfo[i].Name()))
				if err != nil {
					fmt.Println("failed to read directory: " + directory)
					log.Fatal(err)
				}
				for j := 0; j < len(subFileInfo); j++ {
					if filepath2.Ext(subFileInfo[j].Name()) == ext2 {
						_ = os.Remove(filepath2.Join(directory, fileInfo[i].Name(), subFileInfo[j].Name()))
					} else if filepath2.Ext(subFileInfo[j].Name()) == ext1 {
						baseName := strings.Split(subFileInfo[j].Name(), ".")[0]
						convName := baseName + ext2
						basePath := filepath2.Join(directory, fileInfo[i].Name(), subFileInfo[j].Name())
						convPath := filepath2.Join(directory, fileInfo[i].Name(), convName)

						wg.Add(1)
						go obabelWrapper(basePath, convPath, addHydrogens, addCoords, &wg)

					}
				}
			} else if filepath2.Ext(fileInfo[i].Name()) == ext1 {
				baseName := strings.Split(fileInfo[i].Name(), ".")[0]
				convName := baseName + ext2
				basePath := filepath2.Join(directory, fileInfo[i].Name())
				convPath := filepath2.Join(directory, convName)

				wg.Add(1)
				go obabelWrapper(basePath, convPath, addHydrogens, addCoords, &wg)

			} else if filepath2.Ext(fileInfo[i].Name()) == ext2 {
				if deleteOriginal {
					_ = os.Remove(filepath2.Join(directory, fileInfo[i].Name()))
				}
			}
		}
		frame[0] += goRoutinesBatchSize
		frame[1] += goRoutinesBatchSize
		frame[1] = min(frame[1], len(fileInfo)-1)
		wg.Wait()
	}

}

// for file structure directory > file to be converted
func obabelConversion2(directory string, ext1 string, ext2 string, addHydrogens string, addCoords bool) {
	// Read in all files in dir
	fileInfo, err := ioutil.ReadDir(directory)
	if err != nil {
		fmt.Println("failed to read directory: " + directory)
		log.Fatal(err)
	}


	// Iterate through all items in directory
	maxFrame := min(goRoutinesBatchSize,len(fileInfo)-1)
	frame := []int{0,maxFrame}

	for frame[0] < len(fileInfo) {

		wg := sync.WaitGroup{}

		for i := frame[0]; i <= frame[1]; i++ {
			if filepath2.Ext(fileInfo[i].Name()) == ext2 {
				_ = os.Remove(filepath2.Join(directory, fileInfo[i].Name()))
			} else if filepath2.Ext(fileInfo[i].Name()) == ext1 {
				baseName := strings.Split(fileInfo[i].Name(), ".")[0]
				convName := baseName + ext2
				basePath := filepath2.Join(directory, fileInfo[i].Name())
				convPath := filepath2.Join(directory, convName)

				wg.Add(1)
				go obabelWrapper(basePath, convPath, addHydrogens, addCoords, &wg)

			}
		}
		frame[0] += goRoutinesBatchSize
		frame[1] += goRoutinesBatchSize
		frame[1] = min(frame[1], len(fileInfo)-1)
		wg.Wait()
	}
}

func obabelWrapper(path1 string, path2 string, addHydrogens string, addCoords bool, wg *sync.WaitGroup) {
	cmdArgs := []string{path1, "-O", path2}
	if addHydrogens == "add" {
		cmdArgs = append(cmdArgs, "-h")
	} else if addHydrogens == "remove" {
		cmdArgs = append(cmdArgs, "-d")
	}
	if addCoords == true {
		cmdArgs = append(cmdArgs, "--gen3d")
	}

	//fmt.Println(path2)
	//cmdstring := obabel + " -i " + path1 + " -o " + path2
	out, err := exec.Command(obabel, cmdArgs...).CombinedOutput()
	//fmt.Println(string(out))
	wg.Done()
	if err != nil {
		fmt.Println(string(out))
		fmt.Println(err)
		log.Fatal(err)
	}
}

