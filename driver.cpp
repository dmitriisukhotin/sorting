//Name: Dean Sukhotin
//Date: 29/10/2023
//Purpose: Familiarise ourselves with sorting algroithms

#include <iostream>
#include <algorithm>
#include <fstream>
#include <random>
#include <chrono>
#include <iomanip>

using namespace std;

void random(int myData[], int size) {
    srand(time(0));
    for(int i = 0; i < size; i++){
        myData[i] = rand() % 1000000;
    }
}

void insertionSort(int array[], int n, uint64_t &numOfComparisons, uint64_t &numOfSwaps){//main sorting algorithm
    numOfComparisons = 0;
    numOfSwaps = 0;
    for(int i = 1; i < n; i++){
        int key = array[i];
        int j = i - 1;
        while(j >= 0 && array[j] > key){
            numOfComparisons++;
            array[j + 1] = array[j];
            numOfSwaps++;
            j--;
        }
        array[j + 1] = key;
        numOfSwaps++;
    }
}

void selectionSort(int array[], int n, uint64_t &numOfComparisons, uint64_t &numOfSwaps){//main sorting algorithm
    numOfComparisons = 0;
    numOfSwaps = 0;
    for(int i = 0; i < n - 1; i++){
        int minValue = i;
        for(int j = i + 1; j < n; j++){
            numOfComparisons++;
            if (array[j] < array[minValue]){
                minValue = j;
            }
        }
        if(minValue != i){
            numOfSwaps++;
            int temp = array[i];
            array[i] = array[minValue];
            array[minValue] = temp;
        }
    }
}

void bubbleSort(int array[], int n, uint64_t &numOfComparisons, uint64_t &numOfSwaps){//main sorting algorithm
    numOfComparisons = 0;
    numOfSwaps = 0;
    for(int i = 0; i < n - 1; i++){
        for(int j = 0; j < n - i - 1; j++){
            numOfComparisons++;
            if(array[j] > array[j + 1]){
                numOfSwaps++;
                int temp = array[j];
                array[j] = array[j + 1];
                array[j + 1] = temp;
            }
        }
    }
}

void merge(int array[], int left, int middle, int right, uint64_t &numOfComparisons, uint64_t &numOfSwaps){
    numOfComparisons = 0;
    numOfSwaps = 0;
    int newArray1 = middle - left + 1;
    int newArray2 = right - middle;
    int leftArray[newArray1], rightArray[newArray2];
    for(int i = 0; i < newArray1; i++){
        leftArray[i] = array[left + i];
    }
    for(int i = 0; i < newArray2; i++){
        rightArray[i] = array[middle + 1 + i];
    }
    int i = 0, j = 0, k = left;
    while(i < newArray1 && j < newArray2){
        numOfComparisons++;
        if(leftArray[i] <= rightArray[j]){
            array[k] = leftArray[i];
            i++;
        }
        else{
            array[k] = rightArray[j];
            j++;
        }
        k++;
    }
    while(i < newArray1){
        array[k] = leftArray[i];
        i++;
        k++;
    }
    while(j < newArray2){
        array[k] = rightArray[j];
        j++;
        k++;
    }
}

void mergeSort(int array[], int left, int right, uint64_t &numOfComparisons, uint64_t &numOfSwaps){     //extra credit sorting algorithm
    numOfComparisons = 0;
    numOfSwaps = 0;
    if(left < right){
        int middle = (left + (right - left) / 2);
        mergeSort(array, left, middle, numOfComparisons, numOfSwaps);
        mergeSort(array, middle + 1, right, numOfComparisons, numOfSwaps);
        merge(array, left, middle, right, numOfComparisons, numOfSwaps);
    }
}

int partition(int array[], int low, int high, uint64_t &numOfComparisons, uint64_t &numOfSwaps){
    numOfComparisons = 0;
    numOfSwaps = 0;
    int pivot = array[high];
    int i = (low - 1);
    for(int j = low; j <= high - 1; j++){
        numOfComparisons++;
        if(array[j] < pivot){
            i++;
            numOfSwaps++;
            int temp = array[i];
            array[i] = array[j];
            array[j] = temp;
        }
    }
    numOfSwaps++;
    int temp = array[i + 1];
    array[i + 1] = array[high];
    array[high] = temp;
    return (i + 1);
}

void quickSort(int array[], int low, int high, uint64_t &numOfComparisons, uint64_t &numOfSwaps){    //extra credit sorting algorithm
    numOfComparisons = 0;
    numOfSwaps = 0;
    if(low < high){
        int pi = partition(array, low, high, numOfComparisons, numOfSwaps);
        quickSort(array, low, pi - 1, numOfComparisons, numOfSwaps);
        quickSort(array, pi + 1, high, numOfComparisons, numOfSwaps);
    }
}

int main(){
for(int i = 0; i < 10; i++){
    int mySizes[] = {1000, 10000, 100000};
    for(int j = 0; j < 3; j++){
        int mySize = mySizes[j];
        int myData[mySize];
        int sortedData[mySize];
        random(myData, mySize);
        uint64_t numOfComparisons, numOfSwaps;
        double totalCPU1 = 0.0;
        double totalCPU2 = 0.0;
        double totalCPU3 = 0.0;
        double totalCPU4 = 0.0;
        double totalCPU5 = 0.0;
        double totalCPU6 = 0.0;
        double totalCPU7 = 0.0;
        double totalCPU8 = 0.0;
        double totalCPU9 = 0.0;
        double totalCPU10 = 0.0;

        uint64_t totalSwaps1 = 0;
        uint64_t totalSwaps2 = 0;
        uint64_t totalSwaps3 = 0;
        uint64_t totalSwaps4 = 0;
        uint64_t totalSwaps5 = 0;
        uint64_t totalSwaps6 = 0;
        uint64_t totalSwaps7 = 0;
        uint64_t totalSwaps8 = 0;
        uint64_t totalSwaps9 = 0;
        uint64_t totalSwaps10 = 0;

        uint64_t totalComparisons1 = 0;
        uint64_t totalComparisons2= 0;
        uint64_t totalComparisons3 = 0;
        uint64_t totalComparisons4 = 0;
        uint64_t totalComparisons5 = 0;
        uint64_t totalComparisons6 = 0;
        uint64_t totalComparisons7 = 0;
        uint64_t totalComparisons8 = 0;
        uint64_t totalComparisons9 = 0;
        uint64_t totalComparisons10 = 0;

        auto start1 = chrono::high_resolution_clock::now();
        copy(myData, myData + mySize, sortedData);
        bubbleSort(myData, mySize, numOfComparisons, numOfSwaps);

        auto end1 = chrono::high_resolution_clock::now();
        auto time1 = chrono::duration<double>(end1 - start1);
        cout << "CPU time: "<< time1.count() << "s" << endl;
        cout << "Bubble Sort:" << endl;
        cout << "Size: " << mySize << endl;
        cout << "Unsorted data: Number of Comparisons: " << numOfComparisons << ", Number of Swaps: " << numOfSwaps << endl;
        totalCPU1 += time1.count();
        totalSwaps1 += numOfSwaps;
        totalComparisons1 += numOfComparisons;

        auto start2 = chrono::high_resolution_clock::now();
        copy(sortedData, sortedData + mySize, myData);
        insertionSort(myData, mySize, numOfComparisons, numOfSwaps);

        auto end2 = chrono::high_resolution_clock::now();
        auto time2 = chrono::duration<double>(end2 - start2);
        cout << "CPU time: "<< time2.count() << "s" << endl;
        cout << "Insertion Sort:" << endl;
        cout << "Size: " << mySize << endl;
        cout << "Unsorted data: Number of Comparisons: " << numOfComparisons << ", Number of Swaps: " << numOfSwaps << endl;
        totalCPU2 += time2.count();
        totalSwaps2 += numOfSwaps;
        totalComparisons2 += numOfComparisons;

        auto start3 = chrono::high_resolution_clock::now();
        copy(sortedData, sortedData + mySize, myData);
        selectionSort(myData, mySize, numOfComparisons, numOfSwaps);

        auto end3 = chrono::high_resolution_clock::now();
        auto time3 = chrono::duration<double>(end3 - start3);
        cout << "CPU time: "<< time3.count() << "s" << endl;
        cout << "Selection Sort:" << endl;
        cout << "Size: " << mySize << endl;
        cout << "Unsorted data: Number of Comparisons: " << numOfComparisons << ", Number of Swaps: " << numOfSwaps << endl;
        totalCPU3 += time3.count();
        totalSwaps3 += numOfSwaps;
        totalComparisons3 += numOfComparisons;

        auto start4 = chrono::high_resolution_clock::now();
        copy(sortedData, sortedData + mySize, myData);
        mergeSort(myData, 0, mySize - 1, numOfComparisons, numOfSwaps);

        auto end4 = chrono::high_resolution_clock::now();
        auto time4 = chrono::duration<double>(end4 - start4);
        cout << "CPU time: "<< time4.count() << "s" << endl;
        cout << "Merge Sort:" << endl;
        cout << "Size: " << mySize << endl;
        cout << "Unsorted data: Number of Comparisons: " << numOfComparisons << ", Number of Swaps: " << numOfSwaps << endl;
        totalCPU4 += time4.count();
        totalSwaps4 += numOfSwaps;
        totalComparisons4 += numOfComparisons;

        auto start5 = chrono::high_resolution_clock::now();
        copy(sortedData, sortedData + mySize, myData);
        quickSort(myData, 0, mySize - 1, numOfComparisons, numOfSwaps);

        auto end5 = chrono::high_resolution_clock::now();
        auto time5 = chrono::duration<double>(end5 - start5);
        cout << "CPU time: "<< time5.count() << "s" << endl;
        cout << "Quick Sort: " << endl;
        cout << "Size: " << mySize << endl;
        cout << "Unsorted data: Number of Comparisons = " << numOfComparisons << ", Number of Swaps: " << numOfSwaps << endl;
        totalCPU5 += time5.count();
        totalSwaps5 += numOfSwaps;
        totalComparisons5 += numOfComparisons;

        auto start6 = chrono::high_resolution_clock::now();
        copy(sortedData, sortedData + mySize, myData);
        bubbleSort(myData, mySize, numOfComparisons, numOfSwaps);

        auto end6 = chrono::high_resolution_clock::now();
        auto time6 = chrono::duration<double>(end6 - start6);
        cout << "CPU time: "<< time6.count() << "s" << endl;
        cout << "Bubble Sort: " << endl;
        cout << "Size: " << mySize << endl;
        cout << "Sorted data: Number of Comparisons: " << numOfComparisons << ", Number of Swaps: " << numOfSwaps << endl;
        totalCPU6 += time6.count();
        totalSwaps6 += numOfSwaps;
        totalComparisons6 += numOfComparisons;

        auto start7 = chrono::high_resolution_clock::now();
        copy(sortedData, sortedData + mySize, myData);
        insertionSort(myData, mySize, numOfComparisons, numOfSwaps);

        auto end7 = chrono::high_resolution_clock::now();
        auto time7 = chrono::duration<double>(end7 - start7);
        cout << "CPU time: "<< time7.count() << "s" << endl;
        cout << "Insertion Sort: " << endl;
        cout << "Size: " << mySize << endl;
        cout << "Sorted data: Number of Comparisons: " << numOfComparisons << ", Number of Swaps: " << numOfSwaps << endl;
        totalCPU7 += time7.count();
        totalSwaps7 += numOfSwaps;
        totalComparisons7 += numOfComparisons;

        auto start8 = chrono::high_resolution_clock::now();
        copy(sortedData, sortedData + mySize, myData);
        selectionSort(myData, mySize, numOfComparisons, numOfSwaps);

        auto end8 = chrono::high_resolution_clock::now();
        auto time8 = chrono::duration<double>(end8 - start8);
        cout << "CPU time: "<< time8.count() << "s" << endl;
        cout << "Selection Sort: " << endl;
        cout << "Size: " << mySize << endl;
        cout << "Sorted data: Number of Comparisons = " << numOfComparisons << ", Number of Swaps = " << numOfSwaps << endl;
        totalCPU8 += time8.count();
        totalSwaps8 += numOfSwaps;
        totalComparisons8 += numOfComparisons;

        auto start9 = chrono::high_resolution_clock::now();
        copy(sortedData, sortedData + mySize, myData);
        mergeSort(myData, 0, mySize - 1, numOfComparisons, numOfSwaps);

        auto end9 = chrono::high_resolution_clock::now();
        auto time9 = chrono::duration<double>(end9 - start9);
        cout << "CPU time: "<< time9.count() << "s" << endl;
        cout << "Merge Sort: " << endl;
        cout << "Size: " << mySize << endl;
        cout << "Sorted data: Number of Comparisons = " << numOfComparisons << ", Number of Swaps = " << numOfSwaps << endl;
        totalCPU9 += time9.count();
        totalSwaps9 += numOfSwaps;
        totalComparisons9 += numOfComparisons;

        auto start10 = chrono::high_resolution_clock::now();
        copy(sortedData, sortedData + mySize, myData);
        quickSort(myData, 0, mySize - 1, numOfComparisons, numOfSwaps);

        auto end10 = chrono::high_resolution_clock::now();
        auto time10 = chrono::duration<double>(end10 - start10);
        cout << "CPU time: "<< time10.count() << "s" << endl;
        cout << "Quick Sort: " << endl;
        cout << "Size: " << mySize << endl;
        cout << "Sorted data: Number of Comparisons = " << numOfComparisons << ", Number of Swaps = " << numOfSwaps << endl;
        totalCPU10 += time10.count();
        totalSwaps10 += numOfSwaps;
        totalComparisons10 += numOfComparisons;

        double averageCPU1 = totalCPU1 / 10;
        double averageCPU2 = totalCPU2 / 10;
        double averageCPU3 = totalCPU3 / 10;
        double averageCPU4 = totalCPU4 / 10;
        double averageCPU5 = totalCPU5 / 10;
        double averageCPU6 = totalCPU6 / 10;
        double averageCPU7 = totalCPU7 / 10;
        double averageCPU8 = totalCPU8 / 10;
        double averageCPU9 = totalCPU9 / 10;
        double averageCPU10 = totalCPU10 / 10;

        uint64_t averageSwaps1 = totalSwaps1 / 10;
        uint64_t averageSwaps2 = totalSwaps2 / 10;
        uint64_t averageSwaps3 = totalSwaps3 / 10;
        uint64_t averageSwaps4 = totalSwaps4 / 10;
        uint64_t averageSwaps5 = totalSwaps5 / 10;
        uint64_t averageSwaps6 = totalSwaps6 / 10;
        uint64_t averageSwaps7 = totalSwaps7 / 10;
        uint64_t averageSwaps8 = totalSwaps8 / 10;
        uint64_t averageSwaps9 = totalSwaps9 / 10;
        uint64_t averageSwaps10 = totalSwaps10 / 10;

        uint64_t averageComparisons1 = totalComparisons1 / 10;
        uint64_t averageComparisons2 = totalComparisons2 / 10;
        uint64_t averageComparisons3 = totalComparisons3 / 10;
        uint64_t averageComparisons4 = totalComparisons4 / 10;
        uint64_t averageComparisons5 = totalComparisons5 / 10;
        uint64_t averageComparisons6 = totalComparisons6 / 10;
        uint64_t averageComparisons7 = totalComparisons7 / 10;
        uint64_t averageComparisons8 = totalComparisons8 / 10;
        uint64_t averageComparisons9 = totalComparisons9 / 10;
        uint64_t averageComparisons10 = totalComparisons10 / 10;
        
        ofstream myFile;
        if(mySize == 1000){
            myFile.open("1000.csv", ios::app);
            myFile << "Run Time " << "Comparisons " << "Swaps" << endl; 
            myFile << "Unsorted" << endl; 
            myFile << setprecision(11) << averageCPU1 << "," << averageComparisons1 << "," << averageSwaps1 << endl;
            myFile << setprecision(11) << averageCPU2 << "," << averageComparisons2 << "," << averageSwaps2 << endl;
            myFile << setprecision(11) << averageCPU3 << "," << averageComparisons3 << "," << averageSwaps3 << endl;
            myFile << setprecision(11) << averageCPU4 << "," << averageComparisons4 << "," << averageSwaps4 << endl;
            myFile << setprecision(11) << averageCPU5 << "," << averageComparisons5 << "," << averageSwaps5 << endl;
            myFile << "Sorted:" << endl;
            myFile << setprecision(11) << averageCPU6 << "," << averageComparisons6 << "," << averageSwaps6 << endl;
            myFile << setprecision(11) << averageCPU7 << "," << averageComparisons7 << "," << averageSwaps7 << endl;
            myFile << setprecision(11) << averageCPU8 << "," << averageComparisons8 << "," << averageSwaps8 << endl;
            myFile << setprecision(11) << averageCPU9 << "," << averageComparisons9 << "," << averageSwaps9 << endl;
            myFile << setprecision(11) << averageCPU10 << "," << averageComparisons10 << "," << averageSwaps10 << endl;
        }
        else if(mySize == 10000){
            myFile.open("10000.csv", ios::app);
            myFile << "Run Time " << "Comparisons " << "Swaps" << endl; 
            myFile << "Unsorted" << endl; 
            myFile << "Bubble: " << setprecision(11) << averageCPU1 << " , " << averageComparisons1 << " , " << averageSwaps1 << endl;
            myFile << "Insertion: " << setprecision(11) << averageCPU2 << " , " << averageComparisons2 << " , " << averageSwaps2 << endl;
            myFile << "Selection: " << setprecision(11) << averageCPU3 << " , " << averageComparisons3 << " , " << averageSwaps3 << endl;
            myFile << "Merge: " << setprecision(11) << averageCPU4 << " , " << averageComparisons4 << " , " << averageSwaps4 << endl;
            myFile << "Quick: " << setprecision(11) << averageCPU5 << " , " << averageComparisons5 << " , " << averageSwaps5 << endl;
            myFile << "Sorted:" << endl;
            myFile << "Bubble: " << setprecision(11) << averageCPU6 << " , " << averageComparisons6 << " , " << averageSwaps6 << endl;
            myFile << "Insertion: " << setprecision(11) << averageCPU7 << " , " << averageComparisons7 << " , " << averageSwaps7 << endl;
            myFile << "Selection: " << setprecision(11) << averageCPU8 << " , " << averageComparisons8 << " , " << averageSwaps8 << endl;
            myFile << "Merge: " << setprecision(11) << averageCPU9 << " , " << averageComparisons9 << " , " << averageSwaps9 << endl;
            myFile << "Quick: " << setprecision(11) << averageCPU10 << " , " << averageComparisons10 << " , " << averageSwaps10 << endl;
        }
        else if(mySize == 100000){
            myFile.open("100000.csv", ios::app);
            myFile << "Run Time " << "Comparisons " << "Swaps" << endl; 
            myFile << "Unsorted" << endl; 
            myFile << "Bubble: " << setprecision(11) << averageCPU1 << " , " << averageComparisons1 << " , " << averageSwaps1 << endl;
            myFile << "Insertion: " << setprecision(11) << averageCPU2 << " , " << averageComparisons2 << " , " << averageSwaps2 << endl;
            myFile << "Selection: " << setprecision(11) << averageCPU3 << " , " << averageComparisons3 << " , " << averageSwaps3 << endl;
            myFile << "Merge: " << setprecision(11) << averageCPU4 << " , " << averageComparisons4 << " , " << averageSwaps4 << endl;
            myFile << "Quick: " << setprecision(11) << averageCPU5 << " , " << averageComparisons5 << " , " << averageSwaps5 << endl;
            myFile << "Sorted:" << endl;
            myFile << "Bubble: " << setprecision(11) << averageCPU6 << " , " << averageComparisons6 << " , " << averageSwaps6 << endl;
            myFile << "Insertion: " << setprecision(11) << averageCPU7 << " , " << averageComparisons7 << " , " << averageSwaps7 << endl;
            myFile << "Selection: " << setprecision(11) << averageCPU8 << " , " << averageComparisons8 << " , " << averageSwaps8 << endl;
            myFile << "Merge: " << setprecision(11) << averageCPU9 << " , " << averageComparisons9 << " , " << averageSwaps9 << endl;
            myFile << "Quick: " << setprecision(11) << averageCPU10 << " , " << averageComparisons10 << " , " << averageSwaps10 << endl;
            }
            myFile.close();
        }
    }
    return 0;
}
