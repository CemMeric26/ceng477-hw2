#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>

#include "tinyxml2.h"
#include "Triangle.h"
#include "Helpers.h"
#include "Scene.h"
#include <corecrt_math_defines.h>

using namespace tinyxml2;
using namespace std;

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
			{2/(right - left),0,0,-(right + left)/(right - left)},
			{0,2/(top - bottom),0,-(top + bottom)/(top - bottom)},
			{0,0,2/(near - far),-(near + far)/(near - far)},
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
			{0,0,-(far + near)/(far - near),-2*far*near/(far - near)},
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
	// 0 0 1 0.5
	double viewportMatrix[4][4] = {
		{horRes/2,0,0,(horRes-1)/2},
		{0,verRes/2,0,(verRes-1)/2},
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
	double t_e = 0;
	double t_l = 1;
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

	return visible;
}

bool checkBackfaceCulling(Vec4& p1, Vec4& p2, Vec4& p3){

	Vec3 v1 = Vec3(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z); // v1 = p2 - p1
	Vec3 v2 = Vec3(p3.x - p1.x, p3.y - p1.y, p3.z - p1.z); // v2 = p3 - p1

	Vec3 normal = normalizeVec3(crossProductVec3(v1,v2)); // normal = v1 x v2
	Vec3 view = {p1.x,p1.y,p1.z}; // is this true for view vector = view = p1 - origin

	double dotProduct = dotProductVec3(normal, view);

	if(dotProduct > 0){
		return true; // this means that the triangle is not visible
	}
	else{
		return false;
	}

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

			// apply transformations to vertices // first make them Vec4
			Vec4 v1_4 = Vec4(v1->x, v1->y, v1->z, 1,v1->colorId); // first vertex
			Vec4 v2_4 = Vec4(v2->x, v2->y, v2->z, 1, v2->colorId); // second vertex
			Vec4 v3_4 = Vec4(v3->x, v3->y, v3->z, 1, v3->colorId); // third vertex

			// apply transformations to vertices
			v1_4 = multiplyMatrixWithVec4(CameraModelingAndProjectionMatrix, v1_4);
			v2_4 = multiplyMatrixWithVec4(CameraModelingAndProjectionMatrix, v2_4);
			v3_4 = multiplyMatrixWithVec4(CameraModelingAndProjectionMatrix, v3_4);

			// culling will be done here if it is enabled
			if(this->cullingEnabled){
				bool notVisible = checkBackfaceCulling(v1_4, v2_4, v3_4); // returns true if its backfaced
				if(notVisible){
					continue;
				}
			}

			// clipping will be done here
			// NOTE: Clipping will be applied for only Wireframe mode
			// type=0 for wireframe, type=1 for solid
			if(mesh->type == WIREFRAME_MESH){
				// clipping algorithm: liang-barsky
				// returns true if line is visible
				// returns false if line is invisible
				// two edges (v1,v2), (v2,v3), (v3,v1)
				bool visibleLine1 = clippingLiangBarsky(*c1,*c2,v1_4, v2_4, camera->left, camera->right, camera->bottom, camera->top, camera->near, camera->far); 
				bool visibleLine2 = clippingLiangBarsky(*c2,*c3,v2_4, v3_4, camera->left, camera->right, camera->bottom, camera->top, camera->near, camera->far);
				bool visibleLine3 = clippingLiangBarsky(*c3,*c1,v3_4, v1_4, camera->left, camera->right, camera->bottom, camera->top, camera->near, camera->far);

			}
			else{
				// solid mesh
				;
			}

			

			// perspective division before viewport transformation

			// viewport transformation will be done here

			// rasterization will be done here

			//final


	}


}
