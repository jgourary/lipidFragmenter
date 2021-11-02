package main

import "math/rand"

func qsort(a []int, b []string) ([]int, []string) {
	if len(a) < 2 { return a, b }

	left, right := 0, len(a) - 1

	// Pick a pivot
	pivotIndex := rand.Int() % len(a)

	// Move the pivot to the right
	a[pivotIndex], a[right] = a[right], a[pivotIndex]
	b[pivotIndex], b[right] = b[right], b[pivotIndex]

	// Pile elements smaller than the pivot on the left
	for i := range a {
		if a[i] < a[right] {
			a[i], a[left] = a[left], a[i]
			b[i], b[left] = b[left], b[i]
			left++
		}
	}

	// Place the pivot after the last smaller element
	a[left], a[right] = a[right], a[left]
	b[left], b[right] = b[right], b[left]

	// Go down the rabbit hole
	qsort(a[:left], b[:left])
	qsort(a[left + 1:], b[left + 1:])


	return a, b
}

func qsort4(a []string) []string {
	if len(a) < 2 { return a }

	left, right := 0, len(a) - 1

	// Pick a pivot
	pivotIndex := rand.Int() % len(a)

	// Move the pivot to the right
	a[pivotIndex], a[right] = a[right], a[pivotIndex]

	// Pile elements smaller than the pivot on the left
	for i := range a {
		if a[i] < a[right] {
			a[i], a[left] = a[left], a[i]
			left++
		}
	}

	// Place the pivot after the last smaller element
	a[left], a[right] = a[right], a[left]

	// Go down the rabbit hole
	qsort4(a[:left])
	qsort4(a[left + 1:])


	return a
}