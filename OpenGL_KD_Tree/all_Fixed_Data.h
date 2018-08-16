#ifndef _ALL_FIXED_DATA_H
#define _ALL_FIXED_DATA_H


#include "GameObject.h"
#include "Model.h"
#include "GLIncludes.h"

#define number_of_objects 10 // Dont have this number be greater than 1000
#define EPSILON 0.00001

using namespace std;

#pragma region Coord_maker

struct Point_xy
{
	float x;
	float y;

	Point_xy(float u = 0, float v = 0) : x(u), y(v)
	{}
};
struct _Point_xy
{
	float x1;
	float y1;

	_Point_xy(float u = 0.0f, float v = 0.0f) : x1(u), y1(v)
	{}
};

typedef struct Tree_Point
{
	float xpos, ypos;

}Tree_Point;



#pragma endregion

#pragma region program specific Data members

float timestep = 0.016f;
float overlap;
bool collision = false;
int lastState;
const int k = 2;
float min_x = -1.0f;
float min_y = -1.0f;
float max_x = 1.0f;
float max_y = 1.0f;

float min_x_2 = -1.0f;
float min_y_2 = -1.0f;

glm::vec3 pointofcollision;
glm::vec3 minimumTranslationVector;
glm::vec3 mousePos;
glm::vec3 tempVelocity;
glm::vec3 tempPosition;
glm::vec3 tempPosition2;
glm::vec3 tempVelocity2;

float point_body_x[number_of_objects];
float point_body_y[number_of_objects];

float point_tree_x[number_of_objects];
float point_tree_y[number_of_objects];


// vector of scene bodies
std::vector<VertexFormat> lines;
std::vector<GameObject*> bodies;

//vectors for data-trees
std::vector<float> domain;
std::vector<float> range;
std::vector<Tree_Point> Points_Tree;
std::vector<_Point_xy> Point_Coord;

//objects for data-trees
VertexBuffer lineBuffer;
Model* mesh;
Tree_Point point_temp;


#pragma endregion

#endif _ALL_FIXED_DATA_H

