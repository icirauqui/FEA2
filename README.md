# Finite Element Analysis

Program builds a Finite Element Model from point clouds and computes the strain energy.

## Usage

### Compile

<pre class="prettyprint lang-bsh">
<code class="devsite-terminal">git clone https://...</code>
<code class="devsite-terminal">cd ...</code>
<code class="devsite-terminal">mkdir build</code>
<code class="devsite-terminal">cd build</code>
<code class="devsite-terminal">cmake ..</code>
<code class="devsite-terminal">make -j$(nproc)</code>
</pre>


### Run

Parameters:
- path_to_data - relative path to work directory.

<pre class="prettyprint lang-bsh">
<code class="devsite-terminal">./fea2 path_to_data</code>
</pre>


## About

Standalone FEM class from ORB-SLAM2-E at https://github.com/icirauqui/ORB_SLAM2_E.git.
Evaluate results against Abaqus models.

Available element types:
 - C3D6 - Triangular prism
 - C3D8 - Hexadron