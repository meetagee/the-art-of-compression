
/**
 *
 * toqutree (pa3)
 * significant modification of a quadtree .
 * toqutree.cpp
 * This file will be used for grading.
 *
 */

#include "toqutree.h"

toqutree::Node::Node(pair<int,int> ctr, int dim, HSLAPixel a)
	:center(ctr),dimension(dim),avg(a),NW(NULL),NE(NULL),SE(NULL),SW(NULL)
	{}

toqutree::~toqutree(){
	clear(root);
}

toqutree::toqutree(const toqutree & other) {
	root = copy(other.root);
	imgCopy = PNG(other.imgCopy);
}


toqutree & toqutree::operator=(const toqutree & rhs){
	if (this != &rhs) {
		clear(root);
		root = copy(rhs.root);
	}
	return *this;
}

toqutree::toqutree(PNG & imIn, int k){

/* This constructor grabs the 2^k x 2^k sub-image centered */
/* in imIn and uses it to build a quadtree. It may assume  */
/* that imIn is large enough to contain an image of that size. */

/* your code here */
	int width = imIn.width();
	int height = imIn.height();
	cout << "start copy" << endl;
	imgCopy = PNG(imIn);
	if(width != height || (width == height && width != pow(2,k))) {
		PNG* treeImg = new PNG(pow(2,k), pow(2,k));
		int length = pow(2,k);
		for(int i = 0; i < length; i++){
			for(int j= 0; j < length; j++){
				HSLAPixel* p = imIn.getPixel( i + (width - length)/2 , j + (height - length)/2);
				*(treeImg->getPixel(i,j)) = *p;
			}
		}
		root = buildTree(treeImg, k);
		delete treeImg;
		treeImg = NULL;
	}
	else {
		PNG* fullCopy = new PNG(imIn);
		root = buildTree(fullCopy, k);
		delete fullCopy;
		fullCopy = NULL;
	}
}

int toqutree::size() {
/* your code here */
    return count(root);
}


toqutree::Node * toqutree::buildTree(PNG * im, int k) {

/* your code here */

// Note that you will want to practice careful memory use
// In this function. We pass the dynamically allocated image
// via pointer so that it may be released after it is used .
// similarly, at each level of the tree you will want to
// declare a dynamically allocated stats object, and free it
// once you've used it to choose a split point, and calculate
// an average.
	if(k < 0){			//No pixel, return NULL
		return NULL;
	}
    stats* s = new stats(*im);
	//cout << "Clear" << endl;
    pair <int, int> ul(0,0);
    pair <int, int> lr(pow(2,k) - 1, pow(2,k) - 1);
    pair <int, int> center(0, 0);

	if(k == 0){			//reach a leaf, return a node with NULL subnotes
        //delete im;
        //im = NULL;
		Node* ret = new Node(center, k, s->getAvg(ul,lr));
		delete s;
		s = NULL;
		return ret;
	}
	if(k == 1){
		center.first = 1;
		center.second = 1;
		Node* ret = new Node(center, k, s->getAvg(ul,lr));
		//SE
		PNG* imgSE = new PNG(1, 1);
		HSLAPixel * p1 = im->getPixel(1,1);
        *imgSE->getPixel(0,0) = *p1;
		ret->SE = buildTree(imgSE, k-1);
		//NE
		PNG* imgNE = new PNG(1, 1);
		HSLAPixel * p2 = im->getPixel(1,0);
        *imgNE->getPixel(0,0) = *p2;
		ret->NE = buildTree(imgNE, k-1);
		//SW
		PNG* imgSW = new PNG(1, 1);
		HSLAPixel * p3 = im->getPixel(0,1);
        *imgSW->getPixel(0,0) = *p3;
		ret->SW = buildTree(imgSW, k-1);
		//NW
		PNG* imgNW = new PNG(1, 1);
		HSLAPixel * p4 = im->getPixel(0,0);
        *imgNW->getPixel(0,0) = *p4;
		ret->NW = buildTree(imgNW, k-1);
		
		delete imgSE;
		imgSE = NULL;
		delete imgSW;
		imgSW = NULL;
		delete imgNE;
		imgNE = NULL;
		delete imgNW;
		imgNW = NULL;
		delete s;
		s = NULL;
		//free(im);
        //im = NULL;
		return ret;
	} 

	double retEntropy = std::numeric_limits<double>::max();

	int d = pow(2,k);			//length of the matrix
	int a = pow(2,k-1);			//length of subMatrix

	// Finding the optimal spliting point in the green region
    for(int i = pow(2,k-2); i < 3*pow(2,k-2); i++){
        for(int j = pow(2,k-2); j < 3*pow(2,k-2); j++){
			//SE
			pair<int, int> SE_ul(i,j);
			pair<int, int> SE_lr((i+a-1)%d, (j+a-1)%d);
			double entropySE = s->entropy(SE_ul, SE_lr);
			//SW
			int x2 = (i+a)%d;
			int y2 = j;
			pair<int, int> SW_ul(x2, y2);
			pair<int, int> SW_lr((x2 + a - 1)%d, (y2 + a - 1)%d);
			double entropySW = s->entropy(SW_ul, SW_lr);
			//NE
			int x3 = i;
			int y3 = (j+a)%d;
			pair<int, int> NE_ul(x3, y3);
			pair<int, int> NE_lr((x3 + a - 1)%d, (y3 + a - 1)%d);
			double entropyNE = s->entropy(NE_ul, NE_lr);
			// NW
			int x4 = (i+a)%d;
			int y4 = (j+a)%d;
			pair<int, int> NW_ul(x4, y4);
			pair<int, int> NW_lr((x4 + a - 1)%d, (y4 + a - 1)%d);
			double entropyNW = s->entropy(NW_ul, NW_lr);

			double entropy = (entropySE + entropySW + entropyNE + entropyNW)/4;
			if(entropy  < retEntropy) {
				retEntropy = entropy;
				center.first = i;
				center.second = j;
			}
        }
    }
	// Initialize returning node
    Node* ret = new Node(center, k, s->getAvg(ul,lr));

	// Initialize subNodes: NE, NW, SE, SW
	/*_______________SE____________*/
	int SE_x = center.first;
	int SE_y = center.second;
	PNG* imgSE = new PNG(a, a);
    for(int i = 0; i < a; i++){
        for(int j = 0; j < a; j++){
            HSLAPixel * p = im->getPixel((SE_x + i)%d, (SE_y + j)%d);
            *imgSE->getPixel(i,j) = *p;
            p = NULL;           //free
        }
    }
	ret->SE = buildTree(imgSE, k-1);

	/*_______________SW____________*/
	int SW_x = (SE_x + a) % d;
	int SW_y = SE_y;
	PNG* imgSW = new PNG(a, a);
    for(int i = 0; i < a; i++){
        for(int j = 0; j < a; j++){
            HSLAPixel * p = im->getPixel((SW_x + i)%d, (SW_y + j)%d);
            *imgSW->getPixel(i,j) = *p;
             p = NULL;           //free
        }
    }
	ret->SW = buildTree(imgSW, k-1);

	/*_______________NE____________*/
	int NE_x = SE_x;
	int NE_y = (SE_y + a) % d;
	PNG* imgNE = new PNG(a, a);
    for(int i = 0; i < a; i++){
        for(int j = 0; j < a; j++){
            HSLAPixel * p = im->getPixel((NE_x + i)%d, (NE_y + j)%d);
            *imgNE->getPixel(i,j) = *p;
             p = NULL;           //free
        }
    }
	ret->NE = buildTree(imgNE, k-1);

	/*_______________NW____________*/
	int NW_x = (SE_x + a) % d;
	int NW_y = (SE_y + a) % d;
	PNG* imgNW = new PNG(a, a);
    for(int i = 0; i < a; i++){
        for(int j = 0; j < a; j++){
            HSLAPixel * p = im->getPixel((NW_x + i)%d, (NW_y + j)%d);
            *imgNW->getPixel(i,j) = *p;
             p = NULL;           //free
        }
    }
	ret->NW = buildTree(imgNW, k-1);

	/*free the image im after use*/
	delete imgSE;
	imgSE = NULL;
	delete imgSW;
	imgSW = NULL;
	delete imgNE;
	imgNE = NULL;
	delete imgNW;
	imgNW = NULL;
	//free(im);
	//im = NULL;
	
	delete s;
	s = NULL;

	return ret;
}

PNG toqutree::render(){

// My algorithm for this problem included a helper function
// that was analogous to Find in a BST, but it navigated the
// quadtree, instead.

/* your code here */
    //cout << "Hello begin render" << endl;
    
	PNG ret(imgCopy);
	//cout << "HERE" << endl;
	int width = imgCopy.width();
	int height = imgCopy.height();
	int length = pow(2, root->dimension);
	for(int x = 0; x < length; x++){
		for(int y = 0; y < length; y++){
			//cout << "pair: [" << x << "][" << y << "]" << endl;
			HSLAPixel * p = ret.getPixel(x + (width - length)/2, y + (height - length)/2);
			*p = travel(x, y, root);
			p = NULL;
		}
	}
	//PNG ret2 = *ret;
    return ret;
}

HSLAPixel toqutree::travel(int x, int y, Node * node) {
	//cout << "pair: [" << x << "][" << y << "]" << endl;
	//cout << "dimension: " << node->dimension << endl;
    if(node->SE != NULL && node->SW != NULL && node->NE != NULL && node->NW != NULL){           //not the leaf node
        int d = pow(2, node->dimension);
        int a = pow(2, node->dimension - 1);
		//cout << "d: " << d << ", a: " << a << endl;
        //check if the pixel is in SE image
        int SEx1 = node->center.first;
        int SEy1 = node->center.second;
		int SEx2 = (SEx1 + a - 1)%d;
		int SEy2 = (SEy1 + a - 1)%d;
				if(SEx2 >= SEx1 && SEy2 >= SEy1){	// no wrap
					if(SEx1 <= x && x <= SEx2 && SEy1 <= y && y <= SEy2) {
						return travel(x - SEx1, y - SEy1, node->SE);
					}
				}
				else if(SEx2 >= SEx1 && SEy2 < SEy1){ // vertical wrap
					if(SEx1 <= x && SEx2 >= x && 0 <= y && SEy2 >= y) {// upper
                        return travel(x - SEx1, d - SEy1 + y, node->SE);
                    }
                    if(SEx1 <= x && SEx2 >= x && SEy1 <= y && d-1 >= y) { // lower
                        return travel(x - SEx1, y - SEy1, node->SE);
                    }
                }
                else if(SEx2 < SEx1 && SEy2 >= SEy1){	// horizontal wrap
                    if(0 <= x && SEx2 >= x && SEy1 <= y && SEy2 >= y) { // left
                        return travel(d - SEx1 + x, y - SEy1, node->SE);
                    }
                    if(SEx1 <= x && d-1 >= x && SEy1 <= y && SEy2 >= y)	{ // right
                        return travel(x - SEx1, y - SEy1, node->SE);
                    }
                }
                else if(SEx2 < SEx1 && SEy2 < SEy1){ // edge wrap
                    if(0 <= x && x <= SEx2 && 0 <= y && y <= SEy2) { // upper left
                        return travel(d - SEx1 + x, d - SEy1 + y, node->SE);
                    }
                    if(SEx1 <= x && x <= d-1 && 0 <= y && y <= SEy2) {	// upper right
                        return travel(x - SEx1, d - SEy1 + y, node->SE);
                    }
                    if(0 <= x && x <= SEx2 && SEy1 <= y && y <= d-1) {	// lower left
                        return travel(d - SEx1 + x, y - SEy1, node->SE);
                    }
                    if(SEx1 <= x && x <= d-1 && SEy1 <= y && y <= d-1) { // lower right
                        return travel(x - SEx1, y - SEy1, node->SE);
                    }
                }

        //check if the pixel is in SW image
		int SWx1 = (SEx1+a)%d;
		int SWy1 = SEy1;
		int SWx2 = (SWx1 + a - 1)%d;
		int SWy2 = (SWy1 + a - 1)%d;
		//cout << "SWx1: " << SWx1 << ", SWy1: " << SWy1 << ", SWx2: " << SWx2 << ", SWy2: " << SWy2 << endl;
        	if(SWx2 >= SWx1 && SWy2 >= SWy1){	// no wrap
					if(SWx1 <= x && x <= SWx2 && SWy1 <= y && y <= SWy2) {
						return travel(x - SWx1, y - SWy1, node->SW);
					}
				}
				else if(SWx2 >= SWx1 && SWy2 < SWy1){ // vertical wrap
					if(SWx1 <= x && SWx2 >= x && 0 <= y && SWy2 >= y) {// upper
                        return travel(x - SWx1, d - SWy1 + y, node->SW);
                    }
                    if(SWx1 <= x && SWx2 >= x && SWy1 <= y && d-1 >= y) { // lower
                        return travel(x - SWx1, y - SWy1, node->SW);
                    }
                }
                else if(SWx2 < SWx1 && SWy2 >= SWy1){	// horizontal wrap
                    if(0 <= x && SWx2 >= x && SWy1 <= y && SWy2 >= y) { // left
                        return travel(d - SWx1 + x, y - SWy1, node->SW);
                    }
                    if(SWx1 <= x && d-1 >= x && SWy1 <= y && SWy2 >= y)	{ // right
                        return travel(x - SWx1, y - SWy1, node->SW);
                    }
                }
                else if(SWx2 < SWx1 && SWy2 < SWy1){ // edge wrap
                    if(0 <= x && x <= SWx2 && 0 <= y && y <= SWy2) { // upper left
                        return travel(d - SWx1 + x, d - SWy1 + y, node->SW);
                    }
                    if(SWx1 <= x && x <= d-1 && 0 <= y && y <= SWy2) {	// upper right
                        return travel(x - SWx1, d - SWy1 + y, node->SW);
                    }
                    if(0 <= x && x <= SWx2 && SWy1 <= y && y <= d-1) {	// lower left
                        return travel(d - SWx1 + x, y - SWy1, node->SW);
                    }
                    if(SWx1 <= x && x <= d-1 && SWy1 <= y && y <= d-1) { // lower right
                        return travel(x - SWx1, y - SWy1, node->SW);
                    }
                }

        //check if the pixel is in NE image
		int NEx1 = SEx1;
		int NEy1 = (SEy1+a)%d;
		int NEx2 = (NEx1 + a - 1)%d;
		int NEy2 = (NEy1 + a - 1)%d;
		//cout << "NEx1: " << NEx1 << ", NEy1: " << NEy1 << ", NEx2: " << NEx2 << ", NEy2: " << NEy2 << endl;
        	if(NEx2 >= NEx1 && NEy2 >= NEy1){	// no wrap
					if(NEx1 <= x && x <= NEx2 && NEy1 <= y && y <= NEy2) {
						return travel(x - NEx1, y - NEy1, node->NE);
					}
				}
				else if(NEx2 >= NEx1 && NEy2 < NEy1){ // vertical wrap
					if(NEx1 <= x && NEx2 >= x && 0 <= y && NEy2 >= y) {// upper
                        return travel(x - NEx1, d - NEy1 + y, node->NE);
                    }
                    if(NEx1 <= x && NEx2 >= x && NEy1 <= y && d-1 >= y) { // lower
                        return travel(x - NEx1, y - NEy1, node->NE);
                    }
                }
                else if(NEx2 < NEx1 && NEy2 >= NEy1){	// horizontal wrap
                    if(0 <= x && NEx2 >= x && NEy1 <= y && NEy2 >= y) { // left
                        return travel(d - NEx1 + x, y - NEy1, node->NE);
                    }
                    if(NEx1 <= x && d-1 >= x && NEy1 <= y && NEy2 >= y)	{ // right
                        return travel(x - NEx1, y - NEy1, node->NE);
                    }
                }
                else if(NEx2 < NEx1 && NEy2 < NEy1){ // edge wrap
                    if(0 <= x && x <= NEx2 && 0 <= y && y <= NEy2) { // upper left
                        return travel(d - NEx1 + x, d - NEy1 + y, node->NE);
                    }
                    if(NEx1 <= x && x <= d-1 && 0 <= y && y <= NEy2) {	// upper right
                        return travel(x - NEx1, d - NEy1 + y, node->NE);
                    }
                    if(0 <= x && x <= NEx2 && NEy1 <= y && y <= d-1) {	// lower left
                        return travel(d - NEx1 + x, y - NEy1, node->NE);
                    }
                    if(NEx1 <= x && x <= d-1 && NEy1 <= y && y <= d-1) { // lower right
                        return travel(x - NEx1, y - NEy1, node->NE);
                    }
                }

        //check if the pixel is in NW image
		int NWx1 = (SEx1+a)%d;
		int NWy1 = (SEy1+a)%d;
		int NWx2 = (SEx1 + a - 1)%d;
		int NWy2 = (SEy1 + a - 1)%d;
		//cout << "NWx1: " << NWx1 << ", NWy1: " << NWy1 << ", NWx2: " << NWx2 << ", NWy2: " << NWy2 << endl;
        	if(NWx2 >= NWx1 && NWy2 >= NWy1){	// no wrap
					if(NWx1 <= x && x <= NWx2 && NWy1 <= y && y <= NWy2) {
						return travel(x - NWx1, y - NWy1, node->NW);
					}
				}
				else if(NWx2 >= NWx1 && NWy2 < NWy1){ // vertical wrap
					if(NWx1 <= x && NWx2 >= x && 0 <= y && NWy2 >= y) {// upper
                        return travel(x - NWx1, d - NWy1 + y, node->NW);
                    }
                    if(NWx1 <= x && NWx2 >= x && NWy1 <= y && d-1 >= y) { // lower
                        return travel(x - NWx1, y - NWy1, node->NW);
                    }
                }
                else if(NWx2 < NWx1 && NWy2 >= NWy1){	// horizontal wrap
                    if(0 <= x && NWx2 >= x && NWy1 <= y && NWy2 >= y) { // left
                        return travel(d - NWx1 + x, y - NWy1, node->NW);
                    }
                    if(NWx1 <= x && d-1 >= x && NWy1 <= y && NWy2 >= y)	{ // right
                        return travel(x - NWx1, y - NWy1, node->NW);
                    }
                }
                else if(NWx2 < NWx1 && NWy2 < NWy1){ // edge wrap
                    if(0 <= x && x <= NWx2 && 0 <= y && y <= NWy2) { // upper left
                        return travel(d - NWx1 + x, d - NWy1 + y, node->NW);
                    }
                    if(NWx1 <= x && x <= d-1 && 0 <= y && y <= NWy2) {	// upper right
                        return travel(x - NWx1, d - NWy1 + y, node->NW);
                    }
                    if(0 <= x && x <= NWx2 && NWy1 <= y && y <= d-1) {	// lower left
                        return travel(d - NWx1 + x, y - NWy1, node->NW);
                    }
                    if(NWx1 <= x && x <= d-1 && NWy1 <= y && y <= d-1) { // lower right
                        return travel(x - NWx1, y - NWy1, node->NW);
                    }
                }
    }
    else      //reach the leaf
        return node->avg;
    return HSLAPixel();
}

/* oops, i left the implementation of this one in the file! */
void toqutree::prune(double tol){
	pruneHelper(root,tol);
}

void toqutree::pruneHelper(Node* cRoot, double tol){
	if(cRoot->SE != NULL){
        HSLAPixel* p = new HSLAPixel(cRoot->avg.h, cRoot->avg.s, cRoot->avg.l);
        int flag = checkTol(cRoot, tol, p);
        delete p;
        p = NULL;
        if(flag == 0){
            clear(cRoot->SE);
			cRoot->SE = NULL;
            clear(cRoot->SW);
			cRoot->SW = NULL;
            clear(cRoot->SE);
			cRoot->NE = NULL;
            clear(cRoot->NW);
			cRoot->NW = NULL;
        }
		else{
            pruneHelper(cRoot->SE, tol);
            pruneHelper(cRoot->SW, tol);
            pruneHelper(cRoot->NE, tol);
            pruneHelper(cRoot->NW, tol);
        }
	}
}
int toqutree::checkTol(Node* cRoot, double tol, HSLAPixel* p){
    if(cRoot->SE == NULL){  // a node
        double tolerence = p->dist(cRoot->avg);
        if(tolerence >= tol)
            return 1;           //found a leaf with dist > tol, so no prune here
        else
            return 0;
    }
    else{
        int flag = checkTol(cRoot->SE, tol, p);
        if(flag == 1)
            return flag;
		flag = checkTol(cRoot->SW, tol, p);
        if(flag == 1)
            return flag;
		flag = checkTol(cRoot->NE, tol, p);
        if(flag == 1)
            return flag;
		flag = checkTol(cRoot->NW, tol, p);
        if(flag == 1)
            return flag;
        return flag;
    }
}

/* called by destructor and assignment operator*/
void toqutree::clear(Node * & curr){
/* your code here */
    if(curr != NULL){
        if(curr->NE != NULL)
            clear(curr->NE);
        if(curr->NW != NULL)
            clear(curr->NW);
        if(curr->SE != NULL)
            clear(curr->SE);
        if(curr->SW != NULL)
            clear(curr->SW);
        delete curr;
        curr = NULL;
    }
}

/* done */
/* called by assignment operator and copy constructor */
toqutree::Node * toqutree::copy(const Node * other) {
/* your code here */

    if(other != NULL){
        Node * nroot = new Node(other->center, other->dimension, other->avg);
        if(other->NE != NULL)
            nroot->NE = copy(other->NE);
        if(other->NW != NULL)
            nroot->NW = copy(other->NW);
        if(other->SE != NULL)
            nroot->SE = copy(other->SE);
        if(other->SW != NULL)
            nroot->SW = copy(other->SW);
        return nroot;
    }
    else
        return NULL;
}

int toqutree::count(const Node * node) {
    if(node == NULL)
        return 0;
    else
        return 1 + count(node->NE) + count(node->NW) + count(node->SE) + count(node->SW);
}

