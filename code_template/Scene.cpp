#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>

#include "tinyxml2.h"
#include "Triangle.h"
#include "Helpers.h"
#include "Scene.h"

// Include the appropriate headers based on the operating system
#ifdef _WIN32
    #include <direct.h> // For _mkdir on Windows
#else
    #include <sys/stat.h> // For mkdir on Unix-like systems
#endif

using namespace tinyxml2;
using namespace std;

#define M_PI 3.14159265358979323846
/*
	Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *xmlElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *rootNode = xmlDoc.FirstChild();

	// read background color
	xmlElement = rootNode->FirstChildElement("BackgroundColor");
	str = xmlElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	xmlElement = rootNode->FirstChildElement("Culling");
	if (xmlElement != NULL)
	{
		str = xmlElement->GetText();

		if (strcmp(str, "enabled") == 0)
		{
			this->cullingEnabled = true;
		}
		else
		{
			this->cullingEnabled = false;
		}
	}

	// read cameras
	xmlElement = rootNode->FirstChildElement("Cameras");
	XMLElement *camElement = xmlElement->FirstChildElement("Camera");
	XMLElement *camFieldElement;
	while (camElement != NULL)
	{
		Camera *camera = new Camera();

		camElement->QueryIntAttribute("id", &camera->cameraId);

		// read projection type
		str = camElement->Attribute("type");

		if (strcmp(str, "orthographic") == 0)
		{
			camera->projectionType = ORTOGRAPHIC_PROJECTION;
		}
		else
		{
			camera->projectionType = PERSPECTIVE_PROJECTION;
		}

		camFieldElement = camElement->FirstChildElement("Position");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->position.x, &camera->position.y, &camera->position.z);

		camFieldElement = camElement->FirstChildElement("Gaze");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->gaze.x, &camera->gaze.y, &camera->gaze.z);

		camFieldElement = camElement->FirstChildElement("Up");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->v.x, &camera->v.y, &camera->v.z);

		camera->gaze = normalizeVec3(camera->gaze);
		camera->u = crossProductVec3(camera->gaze, camera->v);
		camera->u = normalizeVec3(camera->u);

		camera->w = inverseVec3(camera->gaze);
		camera->v = crossProductVec3(camera->u, camera->gaze);
		camera->v = normalizeVec3(camera->v);

		camFieldElement = camElement->FirstChildElement("ImagePlane");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &camera->left, &camera->right, &camera->bottom, &camera->top,
			   &camera->near, &camera->far, &camera->horRes, &camera->verRes);

		camFieldElement = camElement->FirstChildElement("OutputName");
		str = camFieldElement->GetText();
		camera->outputFilename = string(str);

		this->cameras.push_back(camera);

		camElement = camElement->NextSiblingElement("Camera");
	}

	// read vertices
	xmlElement = rootNode->FirstChildElement("Vertices");
	XMLElement *vertexElement = xmlElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (vertexElement != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = vertexElement->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = vertexElement->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		this->vertices.push_back(vertex);
		this->colorsOfVertices.push_back(color);

		vertexElement = vertexElement->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	xmlElement = rootNode->FirstChildElement("Translations");
	XMLElement *translationElement = xmlElement->FirstChildElement("Translation");
	while (translationElement != NULL)
	{
		Translation *translation = new Translation();

		translationElement->QueryIntAttribute("id", &translation->translationId);

		str = translationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		this->translations.push_back(translation);

		translationElement = translationElement->NextSiblingElement("Translation");
	}

	// read scalings
	xmlElement = rootNode->FirstChildElement("Scalings");
	XMLElement *scalingElement = xmlElement->FirstChildElement("Scaling");
	while (scalingElement != NULL)
	{
		Scaling *scaling = new Scaling();

		scalingElement->QueryIntAttribute("id", &scaling->scalingId);
		str = scalingElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		this->scalings.push_back(scaling);

		scalingElement = scalingElement->NextSiblingElement("Scaling");
	}

	// read rotations
	xmlElement = rootNode->FirstChildElement("Rotations");
	XMLElement *rotationElement = xmlElement->FirstChildElement("Rotation");
	while (rotationElement != NULL)
	{
		Rotation *rotation = new Rotation();

		rotationElement->QueryIntAttribute("id", &rotation->rotationId);
		str = rotationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		this->rotations.push_back(rotation);

		rotationElement = rotationElement->NextSiblingElement("Rotation");
	}

	// read meshes
	xmlElement = rootNode->FirstChildElement("Meshes");

	XMLElement *meshElement = xmlElement->FirstChildElement("Mesh");
	while (meshElement != NULL)
	{
		Mesh *mesh = new Mesh();

		meshElement->QueryIntAttribute("id", &mesh->meshId);

		// read projection type
		str = meshElement->Attribute("type");

		if (strcmp(str, "wireframe") == 0)
		{
			mesh->type = WIREFRAME_MESH;
		}
		else
		{
			mesh->type = SOLID_MESH;
		}

		// read mesh transformations
		XMLElement *meshTransformationsElement = meshElement->FirstChildElement("Transformations");
		XMLElement *meshTransformationElement = meshTransformationsElement->FirstChildElement("Transformation");

		while (meshTransformationElement != NULL)
		{
			char transformationType;
			int transformationId;

			str = meshTransformationElement->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			mesh->transformationTypes.push_back(transformationType);
			mesh->transformationIds.push_back(transformationId);

			meshTransformationElement = meshTransformationElement->NextSiblingElement("Transformation");
		}

		mesh->numberOfTransformations = mesh->transformationIds.size();

		// read mesh faces
		char *row;
		char *cloneStr;
		int v1, v2, v3;
		XMLElement *meshFacesElement = meshElement->FirstChildElement("Faces");
		str = meshFacesElement->GetText();
		cloneStr = strdup(str);

		row = strtok(cloneStr, "\n");
		while (row != NULL)
		{
			int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);

			if (result != EOF)
			{
				mesh->triangles.push_back(Triangle(v1, v2, v3));
			}
			row = strtok(NULL, "\n");
		}
		mesh->numberOfTriangles = mesh->triangles.size();
		this->meshes.push_back(mesh);

		meshElement = meshElement->NextSiblingElement("Mesh");
	}
}

void Scene::assignColorToPixel(int i, int j, Color c)
{
	this->image[i][j].r = c.r;
	this->image[i][j].g = c.g;
	this->image[i][j].b = c.b;
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
	if (this->image.empty())
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<Color> rowOfColors;
			vector<double> rowOfDepths;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(this->backgroundColor);
				rowOfDepths.push_back(1.01);
			}

			this->image.push_back(rowOfColors);
			this->depth.push_back(rowOfDepths);
		}
	}
	else
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				assignColorToPixel(i, j, this->backgroundColor);
				this->depth[i][j] = 1.01;
				this->depth[i][j] = 1.01;
				this->depth[i][j] = 1.01;
			}
		}
	}
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}

/*
	Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
	ofstream fout;

	fout.open(camera->outputFilename.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFilename << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--)
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
*/
void Scene::convertPPMToPNG(string ppmFileName)
{
	string command;

	// TODO: Change implementation if necessary.
	command = "./magick convert " + ppmFileName + " " + ppmFileName + ".png";
	system(command.c_str());
}


// translation transformation is applied in this function
void createTranslationMatrix(Matrix4& T , Translation& translation) {

    T.values[0][3] = translation.tx;
    T.values[1][3] = translation.ty;
    T.values[2][3] = translation.tz;

}
// scaling transformation is applied in this function
void createScalingMatrix(Matrix4& S, Scaling& scaling ){
	S.values[0][0] = scaling.sx;
	S.values[1][1] = scaling.sy;
	S.values[2][2] = scaling.sz;

}

// rotation transformation is applied in this function
void createRotationMatrix(Matrix4& R, Rotation& rotation ){
	double angle = rotation.angle;
	double ux = rotation.ux;
	double uy = rotation.uy;
	double uz = rotation.uz;

	// find orthonormal basis uvw
	// to find v setting the smallest component of u to zero and 
	// swap the other two while negating one:
	Vec3 u,v,w; 

	u = Vec3(ux,uy,uz,0);

	if (ux <= uy && ux <= uz) {
		// ux is the smallest
		v = Vec3(0,-uz,uy,0);
	} else if (uy <= ux && uy <= uz) {
		// uy is the smallest
		v = Vec3(uz,0,-ux,0);		
	} else {
		// uz is the smallest
		v = Vec3(-uy,ux,0,0);
	}

	w = crossProductVec3(u,v);
	v= normalizeVec3(v);
	w = normalizeVec3(w);

	double other[4][4] = { 
		{ux,uy,uz,0},
		{v.x,v.y,v.z,0},
		{w.x, w.y, w.z,0},
		{0,0,0,1}
	};

	Matrix4 M = Matrix4(other);

	double radian = angle * M_PI / 180.0;
	double cosTheta = cos(radian);
	double sinTheta = sin(radian);

	double rotationMatrix[4][4] = {
		{1,0,0,0},
		{0,cosTheta,-sinTheta,0},
		{0,sinTheta,cosTheta,0},
		{0,0,0,1}
	};

	Matrix4 Rx = Matrix4(rotationMatrix);

	double M_inv[4][4] = {
		{ux,v.x,w.x,0},
		{uy,v.y,w.y,0},
		{uz,v.z,w.z,0},
		{0,0,0,1}
	};
	Matrix4 M_inv_matrix = Matrix4(M_inv);
	
	R = multiplyMatrixWithMatrix(M_inv_matrix, multiplyMatrixWithMatrix(Rx, M));

}
// all three transformations are applied in this function
Matrix4 MakeModelingTransformation(Camera *camera, std::vector<Rotation *>& rotations, std::vector<Scaling *>& scalings, std::vector<Translation *>& translations, Mesh& mesh)
{
	Matrix4 M = getIdentityMatrix();
	int numberOfTransformations = mesh.numberOfTransformations; // number of transformations for this mesh


	for(int i=0;i<numberOfTransformations;i++){
		char transformationType = mesh.transformationTypes[i]; // type of transformation
		int transformationId = mesh.transformationIds[i]; // id of transformation

		// first check transformation type and then apply transformation
		// first we apply the last transformation 
		if(transformationType == 'r'){
			Rotation *rotation = rotations[transformationId-1];
			Matrix4 R = getIdentityMatrix();
			
			createRotationMatrix(R, *rotation);

			M = multiplyMatrixWithMatrix(R, M);
		}
		else if(transformationType == 's'){
			Scaling *scaling = scalings[transformationId-1];
			Matrix4 S = getIdentityMatrix();
			
			createScalingMatrix(S, *scaling);

			M = multiplyMatrixWithMatrix(S, M);
		}
		else if(transformationType == 't'){
			Translation *translation = translations[transformationId-1];
			Matrix4 T = getIdentityMatrix(); // init as Identity matrix
			
			createTranslationMatrix(T, *translation);

			M = multiplyMatrixWithMatrix(T, M);

		}

	}

	return M;
}

Matrix4 MakeCameraTransformation(Camera *camera)
{
	Matrix4 M_cam = getIdentityMatrix();

	Vec3 u = camera->u;
	Vec3 v = camera->v;
	Vec3 w = camera->w;

	double cameraMatrix[4][4] = {
		{u.x,u.y,u.z,0},
		{v.x,v.y,v.z,0},
		{w.x,w.y,w.z,0},
		{0,0,0,1}
	};

	Matrix4 C = Matrix4(cameraMatrix);

	double cameraTranslationMatrix[4][4] = {
		{1,0,0,-camera->position.x},
		{0,1,0,-camera->position.y},
		{0,0,1,-camera->position.z},
		{0,0,0,1}
	};

	Matrix4 T = Matrix4(cameraTranslationMatrix);

	M_cam = multiplyMatrixWithMatrix(C, T);

	return M_cam;
}

Matrix4 MakeProjectionTransformation(Camera *camera)
{
	Matrix4 M_proj = getIdentityMatrix();

	double left = camera->left;
	double right = camera->right;
	double bottom = camera->bottom;
	double top = camera->top;
	double near = camera->near;
	double far = camera->far;


	if(camera->projectionType == ORTOGRAPHIC_PROJECTION){ // orthographic projection
		// formula is given in the slides like:
		// 2/(r-l) 0 0 -(r+l)/(r-l)
		// 0 2/(t-b) 0 -(t+b)/(t-b)
		// 0 0 2/(n-f) -(n+f)/(n-f)
		// 0 0 0 1

		double orthographicMatrix[4][4] = {
			{2/(right - left),0,0,(-(right + left))/(right - left)},
			{0,2/(top - bottom),0,(-(top + bottom))/(top - bottom)},
			{0,0,2/(near - far),(-(near + far))/(far - near)},
			{0,0,0,1}
		};

		M_proj = Matrix4(orthographicMatrix);
	}
	else if(camera->projectionType == PERSPECTIVE_PROJECTION){ // perspective projection
		// formula is given in the slides like:
		// 2n/(r-l) 0 (r+l)/(r-l) 0
		// 0 2n/(t-b) (t+b)/(t-b) 0
		// 0 0 -(f+n)/(f-n) -2fn/(f-n)
		// 0 0 -1 0
		double perspectiveMatrix[4][4] = {
			{2*near/(right - left),0,(right + left)/(right - left),0},
			{0,2*near/(top - bottom),(top + bottom)/(top - bottom),0},
			{0,0,(-(far + near))/(far - near),(-2*far*near)/(far - near)},
			{0,0,-1,0}
		};

		M_proj = Matrix4(perspectiveMatrix);
	}

	return M_proj;
}

Matrix4 MakeViewportTransformation(Camera *camera)
{
	Matrix4 M_vp = getIdentityMatrix();

	double horRes = camera->horRes; // horizontal resolution = n_x
	double verRes = camera->verRes; // vertical resolution = n_y

	// matrix formula of Mvp
	// n_x/2 0 0 (n_x-1)/2 + xmin
	// 0 n_y/2 0 (n_y-1)/2 + ymin
	// 0 0 0.5 0.5
	double viewportMatrix[4][4] = {
		{horRes*0.5,0,0,(horRes-1)*0.5},
		{0,verRes*0.5,0,(verRes-1)*0.5},
		{0,0,0.5,0.5},
		{0,0,0,1}
	};

	M_vp = Matrix4(viewportMatrix);

	return M_vp;
}	

bool visibleLB(double den, double num, double& t_e, double& t_l){

	// the formula is given in the slides as bool visible
	double t;
	if(den > 0){ // potentially entering
		t = num/den;
		if(t > t_l){
			return false;
		}
		else if(t > t_e){
			t_e = t; // t_e changed
		}
	}
	else if(den < 0){ // potentially leaving
		t = num/den;
		if(t < t_e){
			return false;
		}
		else if(t < t_l){
			t_l = t; // t_l changed
		}
	}
	else if(num > 0){ // line parallel to edge 
		return false;
	}

	return true;

}

// function for clipping algorithm: liang-barsky choosen
// returns true if line is visible
// returns false if line is invisible
bool clippingLiangBarsky(Color& c1, Color& c2,Vec4& p1, Vec4& p2, double xmin, double xmax, double ymin, double ymax,double zmin, double zmax){
	// You should interpolate the colors for the new vertices. 
	// You need to calculate the new color based on the distances 
	// to the original vertices.
	double t_e = 0.0;
	double t_l = 1.0;
	double dx = p2.x - p1.x;
	double dy = p2.y - p1.y;
	double dz = p2.z - p1.z;
	bool visible = false;
	Color* cComp = new Color();
	cComp->r = c2.r - c1.r; cComp->g = c2.g - c1.g; cComp->b = c2.b - c1.b;


	// code ay be problematic te and t_l values should change in visible function 
	// if so it should be referenced
	// and should we do color calculations here or in the rasterization part
	if(visibleLB(dx,xmin-p1.x,t_e,t_l)){ // left
		if(visibleLB(-dx,p1.x-xmax,t_e,t_l)){ // right
			if(visibleLB(dy,ymin-p1.y,t_e,t_l)){ // bottom
				if(visibleLB(-dy,p1.y-ymax,t_e,t_l)){ // top
					if(visibleLB(dz,zmin-p1.z,t_e,t_l)){ // front
						if(visibleLB(-dz,p1.z-zmax,t_e,t_l)){ // baCK
							if(t_l < 1){
								p2.x = p1.x + t_l*dx;
								p2.y = p1.y + t_l*dy;
								p2.z = p1.z + t_l*dz;
								// color interpolation as in the forum
								c2.r = c1.r + t_l*cComp->r; c2.g = c1.g + t_l*cComp->g; c2.b = c1.b + t_l*cComp->b;
							}
							if(t_e > 0){
								p1.x = p1.x + t_e*dx;
								p1.y = p1.y + t_e*dy;
								p1.z = p1.z + t_e*dz;
								// color interpolation
								c1.r = c1.r + t_e*cComp->r; c1.g = c1.g + t_e*cComp->g; c1.b = c1.b + t_e*cComp->b;
							}
							visible = true;
						}
					}			
				}
			}
		}

	}
	// std::cout << "cComp: " << cComp << std::endl;

	// to avoid memory leak
	delete cComp;

	return visible;
}

bool checkBackfaceCulling(Vec4& p1, Vec4& p2, Vec4& p3){

	Vec3 v1 = Vec3(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z); // v1 = p2 - p1
	Vec3 v2 = Vec3(p3.x - p1.x, p3.y - p1.y, p3.z - p1.z); // v2 = p3 - p1

	Vec3 normal = normalizeVec3(crossProductVec3(v1,v2)); // normal = v1 x v2
	Vec3 view = {p1.x,p1.y,p1.z}; // is this true for view vector = view = p1 - origin
	Vec3 normal_view = normalizeVec3(view); // normalize view vector

	double dotProduct = dotProductVec3(normal, normal_view);

	return dotProduct < 0; // return true if the triangle is facing away from the camera

}


// Helper Functions to update color
void updateColor(Color& color, const Color& colorStep) {
    color.r += colorStep.r;
    color.g += colorStep.g;
    color.b += colorStep.b;
}

// Helper function to calculate color difference.
Color calculateColorDifference(const Color& c1, const Color& c2, double divisor) {
    return { (c2.r - c1.r) / divisor, (c2.g - c1.g) / divisor, (c2.b - c1.b) / divisor };
}

void rasterizeLine(Vec4& p1, Vec4& p2,int y,  Color& initialColor,  Color& dc,
                   std::vector<std::vector<Color>>& image, std::vector<std::vector<double>>& depth, int d, int yIncrement) {
    Color color = initialColor;
    for (int x = p1.x; x <= p2.x;x++) {
        if (depth[x][y] > p1.z) {
            image[x][y] = color;
            depth[x][y] = p1.z;
        }
		if(d * yIncrement < 0){
			y += yIncrement;
			d += (p1.y - p2.y) + (yIncrement * (p2.x - p1.x));
		}
		else{
			d += (p1.y - p2.y);
		}

        updateColor(color, dc);
    
	}
}

void rasterizeLine2(Vec4& p1, Vec4& p2, int x, const Color& initialColor, const Color& dc,
						std::vector<std::vector<Color>>& image, std::vector<std::vector<double>>& depth,int d, int xIncrement) {
    
	Color color = initialColor;
    for (int y = p1.y; y <= p2.y; y++) {
		if (depth[x][y] > p1.z) {
			image[x][y] = color;
			depth[x][y] = p1.z;
		}
		if(d*xIncrement > 0){
			x += xIncrement;
			d += (p2.x - p1.x) + (xIncrement * (p1.y - p2.y));
		}
		else{
			d += (p2.x - p1.x);
		}
		updateColor(color, dc);
	}
}

// Main function line rasterization function
void lineRasterizationFunc(Vec4& p1, Vec4& p2, Color& c1, Color& c2, std::vector<std::vector<Color>>& image, std::vector<std::vector<double>>& depth) {
    double dx = p2.x - p1.x; // x1-x0
    double dy = p2.y - p1.y; // y1-y0

    if (std::abs(dy) <= std::abs(dx)) { // slope is less than 1
        if (p2.x < p1.x) { // if x2 < x1
            std::swap(p1, p2);
            std::swap(c1, c2);
        }
        int yIncrement = (p2.y < p1.y) ? -1 : 1; // if y2 < y1
        int y = p1.y;
		int d = (p1.y - p2.y) + (yIncrement * (p2.x - p1.x)*0.5);

        Color dc = calculateColorDifference(c1, c2, p2.x - p1.x);

        rasterizeLine(p1, p2, y, c1, dc, image, depth, d, yIncrement);
    } else if(std::abs(dy)>std::abs(dx)) {
        if (p2.y < p1.y) { // if y2 < y1
            std::swap(p1, p2);
            std::swap(c1, c2);
        }
        int xIncrement = (p2.x < p1.x) ? -1 : 1; // if x2 < x1
        int x = p1.x; 
		int d = (p2.x - p1.x) + (xIncrement * (p1.y - p2.y)*0.5);
        Color dc = calculateColorDifference(c1, c2, p2.y - p1.y);

        rasterizeLine2(p1, p2, x, c1, dc, image, depth, d, xIncrement);
    }
}


double f_xy(double x, double y, double x0, double y0, double x1, double y1) {
    return x * (y0 - y1) + y * (x1 - x0) + (x0 * y1 - y0 * x1);
}


void triangleRasterizationFunc(Vec4& p0, Vec4& p1, Vec4& p2, Color& c1, Color& c2, Color& c3, std::vector<std::vector<Color> >& image, std::vector<std::vector<double> >& depth, int nx, int ny){
	// triangle rasterization algorithm is given in the slides
	// check according to nx and ny

	int x_min = std::min(std::min(p0.x, p1.x), p2.x);
	x_min = std::max(x_min, 0); // Ensure x_min is not less than 0

	int x_max = std::max(std::max(p0.x, p1.x), p2.x);
	x_max = std::min(x_max, nx - 1); // Ensure x_max is not greater than nx-1

	int y_min = std::min(std::min(p0.y, p1.y), p2.y);
	y_min = std::max(y_min, 0); // Ensure y_min is not less than 0

	int y_max = std::max(std::max(p0.y, p1.y), p2.y);
	y_max = std::min(y_max, ny - 1); // Ensure y_max is not greater than ny-1

	double f01 = f_xy(p2.x, p2.y, p0.x, p0.y, p1.x, p1.y);
	double f12 = f_xy(p0.x, p0.y, p1.x, p1.y, p2.x, p2.y);
	double f20 = f_xy(p1.x, p1.y, p2.x, p2.y, p0.x, p0.y);

	int hor_limit = nx - 1;
	int ver_limit = ny - 1;


	for(int y=y_min; y<= y_max;y++){
		for(int x=x_min; x<= x_max;x++){
			double alpha = f_xy(x,y,p1.x,p1.y,p2.x,p2.y)/f12;
			double beta = f_xy(x,y,p2.x,p2.y,p0.x,p0.y)/f20;
			double gamma = f_xy(x,y,p0.x,p0.y,p1.x,p1.y)/f01;

			if(alpha >= 0 && beta >= 0 && gamma >= 0){
				double z = alpha*p0.z + beta*p1.z + gamma*p2.z;
				if(z < depth[x][y] && y >= 0 && y <= ver_limit  && x >= 0 && x < hor_limit){
					depth[x][y] = z;
					image[x][y].r = round(alpha*c1.r + beta*c2.r + gamma*c3.r);
					image[x][y].g = round(alpha*c1.g + beta*c2.g + gamma*c3.g);
					image[x][y].b = round(alpha*c1.b + beta*c2.b + gamma*c3.b);
				}
			}
		}
			
	}

}


// Helper function to transform vertices
void transformVertices(const Matrix4& AppliedMatrix, Vec3* vertex, Vec4& transformedVertex) {
    // Convert Vec3 to Vec4
    // transformedVertex = Vec4(vertex->x, vertex->y, vertex->z, 1, vertex->colorId);

    // Apply transformation
    transformedVertex = multiplyMatrixWithVec4(AppliedMatrix, transformedVertex);
}

// Function to perform perspective division on a vertex
void perspectiveDivide(double &x, double &y, double &z, double &t) {
    x /= t;
    y /= t;
    z /= t;
    t = 1; // Normalizing the t component after division
}

/*
	Transformations, clipping, culling, rasterization are done here.
*/

void Scene::forwardRenderingPipeline(Camera *camera)
{
	// TODO: Implement this function
	Matrix4 CameraMatrix = MakeCameraTransformation(camera);
	Matrix4 ProjectionMatrix = MakeProjectionTransformation(camera);
	Matrix4 ViewportMatrix = MakeViewportTransformation(camera);

	// now we have to apply all transformations to each vertex/mesh
	int numberOfMeshes = this->meshes.size();

	for(int i=0;i<numberOfMeshes;i++){
		Mesh *mesh = this->meshes[i];
		int numberOfTriangles = mesh->numberOfTriangles;

		// beforehand calculation of the necessary matrices and multiply them to get ONE Transformation matrix
		Matrix4 ModelingMatrix = MakeModelingTransformation(camera, this->rotations, this->scalings, this->translations, *mesh);
		Matrix4 CameraAndModelingMatrix = multiplyMatrixWithMatrix(CameraMatrix, ModelingMatrix);
		Matrix4 CameraModelingAndProjectionMatrix = multiplyMatrixWithMatrix(ProjectionMatrix, CameraAndModelingMatrix);

		for(int j=0;j<numberOfTriangles;j++){
			Triangle triangle = mesh->triangles[j];
			int v1Id = triangle.vertexIds[0];
			int v2Id = triangle.vertexIds[1];
			int v3Id = triangle.vertexIds[2];

			// get vertices and colors
			Vec3 *v1 = this->vertices[v1Id-1];
			Vec3 *v2 = this->vertices[v2Id-1];
			Vec3 *v3 = this->vertices[v3Id-1];

			Color *c1 = this->colorsOfVertices[v1->colorId-1];
			Color *c2 = this->colorsOfVertices[v2->colorId-1];
			Color *c3 = this->colorsOfVertices[v3->colorId-1];

			// Apply transformations to vertices
			Vec4 v1_4 = Vec4(v1->x, v1->y, v1->z, 1, v1->colorId);
			Vec4 v2_4 = Vec4(v2->x, v2->y, v2->z, 1, v2->colorId);
			Vec4 v3_4 = Vec4(v3->x, v3->y, v3->z, 1, v3->colorId);

			transformVertices(CameraModelingAndProjectionMatrix, v1, v1_4);
			transformVertices(CameraModelingAndProjectionMatrix, v2, v2_4);
			transformVertices(CameraModelingAndProjectionMatrix, v3, v3_4);

			// culling will be done here if it is enabled
			if(this->cullingEnabled){
				bool notVisible = checkBackfaceCulling(v1_4, v2_4, v3_4); // returns true if its backfaced
				// std::cout << "notVisible: " << notVisible << std::endl;
				if(notVisible){
					continue;
				}
			}

			// clipping will be done here
			// NOTE: Clipping will be applied for only Wireframe mode
			// type=0 for wireframe, type=1 for solid
			if(mesh->type == WIREFRAME_MESH){
				
				// now perspective division will be done here
				perspectiveDivide(v1_4.x, v1_4.y, v1_4.z, v1_4.t);
				perspectiveDivide(v2_4.x, v2_4.y, v2_4.z, v2_4.t);
				perspectiveDivide(v3_4.x, v3_4.y, v3_4.z, v3_4.t);

				Vec4 v1_4_copy = v1_4;
				Vec4 v2_4_copy = v2_4;
				Vec4 v3_4_copy =v3_4;

				// is it right??
				Color* c1_copy = new Color(); Color* c2_copy = new Color(); Color* c3_copy = new Color();
				c1_copy->r = c1->r; c1_copy->g = c1->g; c1_copy->b = c1->b;
				c2_copy->r = c2->r; c2_copy->g = c2->g; c2_copy->b = c2->b;
				c3_copy->r = c3->r; c3_copy->g = c3->g; c3_copy->b = c3->b;

				Color* c1_copy_2 = new Color(); Color* c2_copy_2 = new Color(); Color* c3_copy_2 = new Color();
				c1_copy_2->r = c1->r; c1_copy_2->g = c1->g; c1_copy_2->b = c1->b;
				c2_copy_2->r = c2->r; c2_copy_2->g = c2->g; c2_copy_2->b = c2->b;
				c3_copy_2->r = c3->r; c3_copy_2->g = c3->g; c3_copy_2->b = c3->b;


				// clipping algorithm: liang-barsky
				// returns true if line is visible
				// returns false if line is invisible
				// two edges (v1,v2), (v2,v3), (v3,v1)
				bool visibleLine1 = clippingLiangBarsky(*c1_copy,*c2_copy, v1_4, v2_4, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
				// std::cout << "visibleLine1: " << visibleLine1 << std::endl;

				bool visibleLine2 = clippingLiangBarsky(*c2_copy_2,*c3_copy, v2_4_copy, v3_4, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
				// std::cout << "visibleLine2: " << visibleLine2 << std::endl;
				bool visibleLine3 = clippingLiangBarsky(*c3_copy_2,*c1_copy_2, v3_4_copy, v1_4_copy, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
				// std::cout << "visibleLine3: " << visibleLine3 << std::endl;

				// now viewport transformation will be done here
				v1_4 = multiplyMatrixWithVec4(ViewportMatrix, v1_4);
				v2_4 = multiplyMatrixWithVec4(ViewportMatrix, v2_4);
				v3_4 = multiplyMatrixWithVec4(ViewportMatrix, v3_4);

				v1_4_copy = multiplyMatrixWithVec4(ViewportMatrix, v1_4_copy);
				v2_4_copy = multiplyMatrixWithVec4(ViewportMatrix, v2_4_copy);
				v3_4_copy = multiplyMatrixWithVec4(ViewportMatrix, v3_4_copy);


				// now rasterization will be done here
				if(visibleLine1){
					// lineRasterization function
					lineRasterizationFunc(v1_4, v2_4, *c1_copy, *c2_copy, this->image, this->depth);
				}
				if(visibleLine2){
					// lineRasterization function
					lineRasterizationFunc(v2_4_copy, v3_4, *c2_copy_2, *c3_copy, this->image, this->depth);
				}
				if(visibleLine3){
					// lineRasterization function
					lineRasterizationFunc(v3_4_copy, v1_4_copy, *c3_copy_2, *c1_copy_2, this->image, this->depth);
				}

				// avoid memory leak
				delete c1_copy; delete c2_copy; delete c3_copy;
				delete c1_copy_2; delete c2_copy_2; delete c3_copy_2;

			}
			else{
				// solid mesh
				// no clipping,Clipping will be applied for only Wireframe mode
				perspectiveDivide(v1_4.x, v1_4.y, v1_4.z, v1_4.t);
				perspectiveDivide(v2_4.x, v2_4.y, v2_4.z, v2_4.t);
				perspectiveDivide(v3_4.x, v3_4.y, v3_4.z, v3_4.t);

				// now viewport transformation will be done here
				v1_4 = multiplyMatrixWithVec4(ViewportMatrix, v1_4);
				v2_4 = multiplyMatrixWithVec4(ViewportMatrix, v2_4);
				v3_4 = multiplyMatrixWithVec4(ViewportMatrix, v3_4);

				Color* c1_copy = new Color(); Color* c2_copy = new Color(); Color* c3_copy = new Color();
				c1_copy->r = c1->r; c1_copy->g = c1->g; c1_copy->b = c1->b;
				c2_copy->r = c2->r; c2_copy->g = c2->g; c2_copy->b = c2->b;
				c3_copy->r = c3->r; c3_copy->g = c3->g; c3_copy->b = c3->b;

				// now rasterization will be done here
				triangleRasterizationFunc(v1_4, v2_4, v3_4, *c1_copy, *c2_copy, *c3_copy, this->image, this->depth, camera->horRes, camera->verRes);
				// avoid memory leak
				delete c1_copy; delete c2_copy; delete c3_copy;


			}

			// avoid memory leak
			// delete v1; delete v2; delete v3;
			// delete c1; delete c2; delete c3;

		}

	}

}
