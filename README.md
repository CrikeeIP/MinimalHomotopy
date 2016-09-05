# MinimalHomotopy
["Measuring Similarity Between Curves on 2-Manifolds via Homotopy Area"] (https://github.com/CrikeeIP/MinimalHomotopy/blob/master/background/ChambersWang_article.pdf) is a paper published by Erin Wolf Chambers and Yusu Wang in 2013. Here we apply their results to polygons in the plane to measure their similarity.

##Introduction
This repository is home to a **C++** implementation of the algorithm described by Chambers and Wang.
The aim is to measure the similarity between two given polygons **P**, **Q**  in the plane R^2. This measure will be the *minimal homotopy area*, i.e. the minimal area to be swept while continuously deforming **P** into **Q** (or vice versa). Thus, two polygons are similar to each other with respect to this measure, if the area needed to deform one into the other is small. 

To that end, we investigate the closed curve **l := P°rev(Q)** (where '**rev(Q)**' is is simply **Q** reversed, and '**°**' is the concatenation). The curve **l** divides the plane into several areas or "faces" **f_i**. Let **A(f_i)** be the area of **f_i** and **w(f_i)** be the winding number of **l** around face **f_i**. As Chambers/Wang proved, the sum of the products **A(f_i)*w(f_i)** can be minimized in a certain way to obtain the minimal homotopy area.  
Given two polygons fulfilling the restrictions, the algorithm performs this minimization and outputs the minimal homotopy area.


##Restrictions
1. Both polygons **P**, **Q** must be non self-intersecting and
2. have to share the same starting and ending points.
3. All intersections between **P** and **Q** must me transversal


##Dependencies
Two lightweight header-only libraries:  
1. [Geometry](https://github.com/CrikeeIP/Geometry)  
2. [FunctionalPlus](https://github.com/Dobiasd/FunctionalPlus)  

And boost (::graph, to be exact)  
4. [Boost](http://www.boost.org/)

##Installation
Before the installation, make sure you have installed [Boost](http://www.boost.org/).
Subsequently, you can use one of the following two alternative ways to install MinimalHomotopy:

***Alternative 1:***  
Clone this repository
```sh
git clone https://github.com/CrikeeIP/MinimalHomotopy
```
and execute the the *install.sh* script delivered with it.

***Alternative 2:***  
If you're uncomfortable running who-knows-what-they'll-do foreign scripts (and are too tired to check them before execution), you can also do it manually:
```sh
#Clone the reporitory
git clone https://github.com/CrikeeIP/OPTICS-Clustering
#Download CImg header
cd OPTICS-Clustering
cd include
cd optics
mkdir CImg
cd Cimg
wget https://raw.githubusercontent.com/dtschump/CImg/master/CImg.h
cd ..
#Download Geometry header
mkdir Geometry
cd Geometry
wget https://raw.githubusercontent.com/CrikeeIP/Geometry/master/include/geometry/geometry.h
cd ..
cd ..
cd ..

#Clone the FunctionalPlus repository and install it
git clone https://github.com/Dobiasd/FunctionalPlus
cd FunctionalPlus
mkdir build
cd build
cmake ..
sudo make install
cd ..
cd ..

#Run the test
cd OPTICS-Clustering
cd test
g++ --std=c++11 -I../include main.cpp -lX11 -lpthread
./a.out
```



##License
Distributed under the *MIT Software License* (X11 license). (See accompanying file LICENSE.)
