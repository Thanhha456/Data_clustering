Project for assessment in "Algorithmic Thinking (Part 2)" course at Rice university.   
https://www.coursera.org/learn/algorithmic-thinking-2  
## Assessment: *Clustering Algorithms and Comparision of Clustering Algorithms*
#Target:  Implementing two methods for computing closest pairs and two methods for clustering data and then comparing these two clustering methods in terms of efficiency, automation, and quality

Input:
- Each data set corresponds to a county in the USA (identified by a unique 5 digit string called a FIPS county code) and includes information on the total population of the county and its per-capita lifetime cancer risk due to air toxics.
The raw data set includes 3108 counties. Taking the thresholds 3.5, 4.5, and 5.5 yields smaller data sets with 896, 290, and 111 counties, respectively. These four data sets will be our primary test data for your clustering methods and are available for download here in CSV form: ( 3108 counties, 896 counties, 290 counties, 111 counties).
- Provided  Cluster(FIPS_codes, horiz_center, vert_center, total_population, average_risk) and two methods distance(other_cluster) and merge_clusters(other_cluster) 

Outcomes:
* Efficiency *
 - The running times of the functions slow_closest_pair and fast_closest_pair should be O(n^2) and O(nlog^2(n)), respectively, which should be easily noticed by the shape of their graphics. The running time of fast closest pairs is almost linear, while slow closest pairs has quadratic runnning time.
- Assuming that hierarchical_clustering uses fast closest pair and that k-means clustering always uses a small fixed of iterations we can say that k-means clustering definetely faster than hierarchical clustering. As fast_closest_pair has the running times O(nlog^2(n)) then hierarchical_clustering would run O(n^2log^2(n)). For k-means_clustering the running time is either O(n) or O(n^2) depending on whether the size of the  output clusters is fixed or varries on the input size.

* Automation and Quality *
(Which method requires less human supervision with less error)
- Base on the illustrations of two methods of clustering, we can see that not all clusterings are the same, some clusterings are better than others. Given a list of clusters, the distortion of the clustering is the sum of the errors associated with its clusters.
 From the graphics of the distortions on difference numbers of the output size of both methods of clustering, we can see that on 111 counties, smaller data, hierarchical clustering has lower distortion,and  with the rise of the input size, k-means clustering has clearly smaller distortion than hierarchical clustering, so no method is dominated on all areas.
- From the other side, k-means clustering depends on the conception of choising the initial centers, which leads to some weeknesses of the method. Let's have a look at the illustrations of the methods on the same data 3108 counties of the US cancer risk at the West coast side.The hierarchical clustering has different shapes and centers than the othe one. The k-means clustering has higher distortion because choosed four centers, one in Washington state, but three in Caliphornia, two of them from south California. The reason is the larger populations in these three counties. So in term of quality on large data, no method is determined.
