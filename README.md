# Spectral_K-means
In this project I implemented a version of the normalized spectral clustering algorithm.
The algorithm's purpose is to find k clusters that extract data with different amound of features into k vectors only. 

This project is based on 2 previous smaller projects:
1. Implementation of the popular K-means algorithm (in both C & Python indepndently)
2. Implementation of K-means ++ algorithm (an upgrade for the previous algorithm - based on find k initial centroids and than start performing K-means algorithm)

In this project I handled writing code in 2 different languages - C & Python.
I also wrote an extension - Python C API.
This extension links between both languages and serves as an arithmetic & logic improvment for the project.
The extension is used for sending data from Python to C, perform calculations based on C code and finally return the output back to Python.

This project performs different procedures (independently) from the noramlized spectral clustering algorithm.
As a user, you can call a specific procedure both from Python and C.

You can also diagonalize a real-symmetric matrix via the Jacobi procedure for diagonalization (as a stand-alone procedure).

This project emphasizes the importance of working flow between 2 different coding languages.
In Python - I used the libraries: Pandas & NumPy.
Furthermore, I used a module written in C for performing most of the arithmetic calculations vis C instead of Python (time improvment).
