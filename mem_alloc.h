#include <iostream>

// allocation for 1D arrays
template <class Type> 
void allocate(Type *&array, int size1)
{
    array = new Type[size1];
    
    if (array == NULL) {
        std::cout << "Failed to allocate memory!" << std::endl;
        exit(1); 
    }
    
    for (int i = 0; i < size1; i++) {
        array[i] = 0.0;
    }
}

// memory allocation for 2D arrays
template <class Type>
void allocate(Type **&array, int size1, int size2)
{
    array = new Type*[size1];
    if (array == NULL) {
        std::cout << "Failed to allocate memory!" << std::endl;
        exit(1);
    }
    
    // allocate each row and initialize to 0
    for (int i = 0; i < size1; i++) {
        array[i] = new Type[size2];
        if (array[i] == NULL) {
            std::cout << "Failed to allocate memory!" << std::endl;
            exit(1);
        }
        
        for (int j = 0; j < size2; j++) {
            array[i][j] = 0.0;
        }
    }
}

// deallocation for 1D arrays
template <class Type>
void deallocate(Type *array, int size1)
{
    delete[] array;
}

// deallocation for 2D arrays
template <class Type>
void deallocate(Type **array, int size1, int size2)
{
    for (int i = 0; i < size1; i++) {
        delete[] array[i];
    }
    delete[] array;
}
