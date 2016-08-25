# MinimalHomotopy
["Measuring Similarity Between Curves on 2-Manifolds via Homotopy Area"] (https://github.com/CrikeeIP/MinimalHomotopy/blob/master/background/ChambersWang_article.pdf) is a paper published by Erin Wolf Chambers and Yusu Wang in 2013. Here we apply their results to polygons in the plane.

##Introduction
This repository is home to a **C++** implementation of the algorithm as described by Chambers and Wang.
The aim is, to measure the similarity between two curves (i.e. polygons) **P**, **Q**  in the plane.

To that end, we investigate the closed curve **P°rev(Q)** (where rev(Q) is is simply Q reversed).


Restrictions are, that both curves be non self-intersecting and share the same start and end point.

For further explanation on how the algorithm works, see the original [paper] (https://github.com/CrikeeIP/MinimalHomotopy/blob/master/background/ChambersWang_article.pdf). 

##Usage
Suppose you have a set of points in R^n, described in cartesion coordinates, and wonder if they have a cluster structure.
Then you might consider using this library, as it offers an interface that lets you draw a [reachability-plot](https://github.com/CrikeeIP/OPTICS-Clustering/blob/master/resources/reachabilityplot.png) with two lines of code:

```cpp
#include <optics.h>

typedef std::vector<double> point; //A list of n cartesian coordinates makes a point
std::vector<point> points; //Your list of points goes here

int main(){
   auto reach_dists = optics::compute_reachability_dists( points, 10, 100 );
   optics::draw_reachability_plot( reach_dists, "D:/reachdists.bmp" );
}
```

##Restrictions

1. Both curves must be non self-intersecting and
2. have to share the same start and end point.

##Dependencies
Three lightweight header-only libraries:  
1. [Geometry](https://github.com/CrikeeIP/Geometry)  
2. [FunctionalPlus](https://github.com/Dobiasd/FunctionalPlus)  
3. [CImg](https://github.com/dtschump/CImg)

And boost (::geometry and ::index, to be exact)  
4. [Boost](http://www.boost.org/)


##Disclaimer

The functionality in this library initially grew due to my personal need for it while using C++ on a regular basis. I try my best to make it error free and as comfortable to use as I can. The API still might change in the future. If you have any suggestions, find errors, miss some functions or want to give general feedback/criticism, I'd love to hear from you. Of course, [contributions](https://github.com/CrikeeIP/OPTICS-Clustering/pulls) are also very welcome.

##License

Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)
