# MinimalHomotopy
["Measuring Similarity Between Curves on 2-Manifolds via Homotopy Area"] (https://github.com/CrikeeIP/MinimalHomotopy/blob/master/background/ChambersWang_article.pdf) is a paper published by Erin Wolf Chambers and Yusu Wang in 2013. Here we apply their results to polygons in the plane to measure their similarity.

##Introduction
This repository is home to a **C++** implementation of the algorithm described by Chambers and Wang.
The aim is to measure the similarity between two given polygons **P**, **Q**  in the plane R^2. This measure will be the *minimal homotopy area*, i.e. the minimal area to be swept while continuously deforming **P** into **Q** (or vice versa). Thus, two polygons are similar to each other with respect to this measure, if the area needed to deform one into the other is small. 

To that end, we investigate the closed curve **l := P°rev(Q)** (where '**rev(Q)**' is is simply **Q** reversed, and '**°**' is the concatenation). The curve **l** divides the plane into several areas or "faces" **f_i**. Let **A(f_i)** the area of **f_i** and **w(f_i)** the winding number of **l** around face **f_i**. As Chambers/Wang proved, the sum over the products **A(f_i)*w(f_i)** can be minimized in a certain way to obtain the minimal homotopy area.


##Restrictions
1. Both polygons **P**, **Q** must be non self-intersecting and
2. have to share the same start and end point.
3. All intersections between **P** and **Q** must me trensversal


##Dependencies
Two lightweight header-only libraries:  
1. [Geometry](https://github.com/CrikeeIP/Geometry)  
2. [FunctionalPlus](https://github.com/Dobiasd/FunctionalPlus)  

And boost (::geometry and ::index, to be exact)  
4. [Boost](http://www.boost.org/)


##Disclaimer
The functionality in this library initially grew due to my personal need for it while using C++ on a regular basis. I try my best to make it error free and as comfortable to use as I can. The API still might change in the future. If you have any suggestions, find errors, miss some functions or want to give general feedback/criticism, I'd love to hear from you. Of course, [contributions](https://github.com/CrikeeIP/OPTICS-Clustering/pulls) are also very welcome.

##License
Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)
