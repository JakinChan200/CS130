#include <algorithm>
#include "hierarchy.h"
#include <queue>

// Reorder the entries vector so that adjacent entries tend to be nearby.
// You may want to implement box.cpp first.
void Hierarchy::Reorder_Entries()
{
    if(!entries.size()) return;
    //if(entries.size() <= 8) return;

    //Split the objects into 4 groups based on x value
    unsigned int pivotx1 = Partition(0, 0, entries.size()-1);
    unsigned int pivotx2 = Partition(0, 0, pivotx1);
    unsigned int pivotx3 = Partition(0, pivotx2+1, entries.size()-1);

    //Split the objects into 8 goups based on y value
    unsigned int pivoty1 = Partition(1, 0, pivotx2);
    unsigned int pivoty2 = Partition(1, pivotx2+1, pivotx1);
    unsigned int pivoty3 = Partition(1, pivotx1+1, pivotx3);
    unsigned int pivoty4 = Partition(1, pivotx3+1, entries.size()-1);

    //Split the objects into 16 goups based on y value
    Partition(2, 0, pivoty1);
    Partition(2, pivoty1+1, pivotx2);
    Partition(2, pivotx2+1, pivoty2);
    Partition(2, pivoty2+1, pivotx1);
    Partition(2, pivotx1+1, pivoty3);
    Partition(2, pivoty3+1, pivotx3);
    Partition(2, pivotx3+1, pivoty4);
    Partition(2, pivoty4+1, entries.size()-1);
}

// Populate tree from entries.
void Hierarchy::Build_Tree()
{
    if(!entries.size()) return;
    tree.resize((2*entries.size())-1); //make sure tree is the correct size

    //populate the tree
    int difference = tree.size()-entries.size();
    for(int i = tree.size()-1; i >= 0; i--){
        if(i >= difference){
            tree[i] = entries[i - difference].box;
        }

        if(i % 2 == 1){
            if((long unsigned int)i != tree.size()-1){
                tree[(i-1)/2] = tree[i];
            }else{
                tree[(i-1)/2] = tree[i].Union(tree[i+1]);
            }
        }
    }
}

// Return a list of candidates (indices into the entries list) whose
// bounding boxes intersect the ray.
void Hierarchy::Intersection_Candidates(const Ray& ray, std::vector<int>& candidates) const
{
    int difference = tree.size() - entries.size();
    std::queue<int> intersect;

    if(tree[0].Intersection(ray)){
        intersect.push(0);
    }

    int cur;
    int child;
    while(!intersect.empty()){
        cur = intersect.front();
        intersect.pop();
        if(cur < difference){
            child = (cur*2)+1;
            if((long unsigned int)child < tree.size() && tree[child].Intersection(ray))
                intersect.push(child);
            child++;
            if((long unsigned int)child < tree.size() && tree[child].Intersection(ray))
                intersect.push(child);
        }else{
            candidates.push_back(cur - difference);
        }
    }
}

int Hierarchy::Partition(int axis, int lower, int higher){
    if(higher - lower <= 1){ //lower >= higher?
        return lower;
    }

    //int pivot = average of 3 numbers;
    int pivot = (entries[lower].box.lo[axis] + entries[higher].box.lo[axis] + entries[(higher+lower)/2].box.lo[axis])/3;

    int lowerIndex = lower-1;
    for(int i = lower+1; i <= higher; i++){
        if(entries[i].box.lo[axis] < pivot){
            lowerIndex++;
            std::swap(entries[lowerIndex], entries[i]);
            // Entry temp = entries[lowerIndex];
            // entries[lowerIndex] = entries[i];
            // entries[i] = temp;
        }
    }
    return lowerIndex;
}