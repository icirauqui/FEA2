# Finite Element Analysis

The repository contains a self-contained class "FEA" for performing Finite Element Analysis. Currently the library supports the following 3D elements:
- Triangular Prism: C3D6.
- Hexadron: C3D8.

The program assembles the Stiffness Matrix and computes the Strain Energy using the known Force vector. 

The only dependency is Eigen3, used for linear algebra operations.

The goal of the class is to be integrated in other projects, with this repository being used for testing purposes only. Information about elements and connectivity must be provided beforehand, we present examples of the data required for each element in the data folder.


## Usage

### Compile

<pre class="prettyprint lang-bsh">
<code class="devsite-terminal">git clone https://github.com/icirauqui/FEA2.git</code>
<code class="devsite-terminal">cd FEA2</code>
<code class="devsite-terminal">mkdir build</code>
<code class="devsite-terminal">cd build</code>
<code class="devsite-terminal">cmake ..</code>
<code class="devsite-terminal">make -j$(nproc)</code>
</pre>


### Run

<pre class="prettyprint lang-bsh">
<code class="devsite-terminal">./fea2</code>
</pre>


## About

The class has been used for processing non-rigid environments in Computer Vision algorithms:
 - Visual SLAM: [ORB-SLAM2-E](https://github.com/icirauqui/ORB_SLAM2_E.git).
 - Alignment of Point Clouds from NRSfM: [PCC](https://github.com/icirauqui/pc_compare.git).

In the former we present a validation of the class results against Abaqus, a commercial Finite Element Analysis software. 
