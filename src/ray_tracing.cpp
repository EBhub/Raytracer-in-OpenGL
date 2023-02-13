#include "ray_tracing.h"
#include "disable_all_warnings.h"
// Suppress warnings in third-party code.
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/vector_relational.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>
#include <iostream>
#include <limits>

bool pointInTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& n, const glm::vec3& p) {
	// creates plane with triangle vertices
	Plane plane = trianglePlane(v0, v1, v2);

	// gets sides of triangles
	glm::vec3 v1v0 = v1 - v0;
	glm::vec3 v2v1 = v2 - v1;
	glm::vec3 v0v2 = v0 - v2;

	// gets vector of each vertex to the point
	glm::vec3 A0 = p - v0;
	glm::vec3 A1 = p - v1;
	glm::vec3 A2 = p - v2;

	// gets dotproduct of each perpendicular vector to the triangles normal
	float a = glm::dot(n, glm::cross(v1v0, A0));
	float b = glm::dot(n, glm::cross(v2v1, A1));
	float c = glm::dot(n, glm::cross(v0v2, A2));

	// checks if the point is to the left of (or on) all the edges
	return (a >= 0 && b >= 0 && c >= 0);
}

bool intersectRayWithPlane(const Plane& plane, Ray& ray) {
	// get distance normal to origin
	float dot = glm::dot(ray.direction, plane.normal);
	// calculate new t from ray
	float t = (plane.D - glm::dot(ray.origin, plane.normal)) / dot;

	// check is dotproduct is positive, bc it cannot be negative or else different facing direction
	// check if t > 0 else same issue
	// check if new t is less than old t = new point is closer
	if (abs(dot) >= 1e-6 && t >= 0 && ray.t > t) {
		// assign new t to ray
		ray.t = t;
		return true;
	}
	else {
		return false;
	}
}

Plane trianglePlane(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2) {
	Plane plane;
	// get edges of triangle (only 2 needed for dotproduct so no need for third edge)
	glm::vec3 v1v0 = (v1 - v0);
	glm::vec3 v2v0 = (v2 - v0);

	// get normal of triangle aka plane
	glm::vec3 N = glm::cross(v1v0, v2v0);
	// normalize bc everything has to be normalized :)
	N = glm::normalize(N);

	// get distance normal to origin
	float D = glm::dot(N, v0);

	// assign plane's new parameters
	plane.D = D;
	plane.normal = N;

	return plane;
}

/// Input: the three vertices of the triangle
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, Ray& ray, HitInfo& hitInfo) {
	// creates plane with triangle's vertices
	Plane plane = trianglePlane(v0, v1, v2);
	// get ray.t for rollback
	float tempT = ray.t;
	bool intsecRayPlane = intersectRayWithPlane(plane, ray);

	// check if ray intersects plane
	if (intsecRayPlane) {
		// get intersection point
		glm::vec3 P = ray.origin + ray.t * ray.direction;
		// check if point is inside triangle
		bool isPointInTriangle = pointInTriangle(v0, v1, v2, plane.normal, P);

		if (isPointInTriangle) {
			hitInfo.normal = plane.normal;
			return true;
		}
	}

	// rollback ray.t
	ray.t = tempT;

	return false;
}

/// Input: a sphere with the following attributes: sphere.radius, sphere.center
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const Sphere& sphere, Ray& ray, HitInfo& hitInfo) {
	//ray1.direction = glm::normalize(ray1.direction);
	glm::vec3 location = ray.origin - sphere.center;

	//directon is normalized thus a = 1, (the length)
	float a = glm::length(ray.direction);
	float b = 2.0f * glm::dot(ray.direction, location);
	float c = glm::dot(location, location) - pow(sphere.radius, 2);

	float determinant = glm::pow(b, 2) - 4 * a * c;

	if (determinant < 0) {
		return false;
	}
	if (determinant == 0) {
		ray.t = -(b / (2 * a));
		return true;
	}

	float t0 = (-b + sqrt(determinant)) / (2 * a);
	float t1 = (-b - sqrt(determinant)) / (2 * a);

	if (t0 < 0 && t1 < 0) {
		return false;
	}

	if (t0 < 0) {
		t0 = std::numeric_limits<float>::max();
	}
	if (t1 < 0) {
		t1 = std::numeric_limits<float>::max();
	}

	if (t0 > t1) {
		ray.t = t1;
	}
	else {
		ray.t = t0;
	}

	// get normal
	glm::vec3 P = ray.origin + ray.direction * ray.t;
	glm::vec3 N = P - sphere.center;
	N = glm::normalize(N);

	hitInfo.normal = N;
	hitInfo.material = sphere.material;
}
/// Input: a sphere with the following attributes: sphere.radius, sphere.center
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const SphericalLight& sphere, Ray& ray, HitInfo& hitInfo) {
	//ray1.direction = glm::normalize(ray1.direction);
	glm::vec3 location = ray.origin - sphere.position;

	//directon is normalized thus a = 1, (the length)
	float a = glm::length(ray.direction);
	float b = 2.0f * glm::dot(ray.direction, location);
	float c = glm::dot(location, location) - pow(sphere.radius, 2);

	float determinant = glm::pow(b, 2) - 4 * a * c;

	if (determinant < 0) {
		return false;
	}
	if (determinant == 0) {
		ray.t = -(b / (2 * a));
		return true;
	}

	float t0 = (-b + sqrt(determinant)) / (2 * a);
	float t1 = (-b - sqrt(determinant)) / (2 * a);

	if (t0 < 0 && t1 < 0) {
		return false;
	}

	if (t0 < 0) {
		t0 = std::numeric_limits<float>::max();
	}
	if (t1 < 0) {
		t1 = std::numeric_limits<float>::max();
	}

	if (t0 > t1) {
		ray.t = t1;
	}
	else {
		ray.t = t0;
	}
}

/// Input: an axis-aligned bounding box with the following parameters: minimum coordinates box.lower and maximum coordinates box.upper
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const AxisAlignedBox& box, Ray& ray) {
	// X ---------------------------------------------------------------
	// get t_min/t_max for x-axis
	float t_xmin = (box.lower.x - ray.origin.x) / ray.direction.x;
	float t_xmax = (box.upper.x - ray.origin.x) / ray.direction.x;

	// get t_in and t_out for x-axis
	float t_inx = std::min(t_xmax, t_xmin);
	float t_outx = std::max(t_xmax, t_xmin);

	// Y ---------------------------------------------------------------
	// get t_min/t_max for y-axis
	float t_ymin = (box.lower.y - ray.origin.y) / ray.direction.y;
	float t_ymax = (box.upper.y - ray.origin.y) / ray.direction.y;

	// get t_in and t_out for y-axis
	float t_iny = std::min(t_ymax, t_ymin);
	float t_outy = std::max(t_ymax, t_ymin);

	// check if t_inx > t_outy
	// OR check if t_iny > t_outx
	// return false, bc missed hit
	if ((t_inx > t_outy) || (t_iny > t_outx)) {
		return false;
	}

	// checks if t_iny > t_inx
	// if true, t_inx = t_iny, bc new lowest t_in
	if (t_iny > t_inx) {
		t_inx = t_iny;
	}

	// checks if t_outy > t_outx
	// if true, t_outx = t_outy, bc new highest t_out
	if (t_outx > t_outy) {
		t_outx = t_outy;
	}

	// Z ---------------------------------------------------------------
	// get t_min/t_max for z-axis
	float t_zmin = (box.lower.z - ray.origin.z) / ray.direction.z;
	float t_zmax = (box.upper.z - ray.origin.z) / ray.direction.z;

	// get t_in and t_out for y-axis
	float t_inz = std::min(t_zmax, t_zmin);
	float t_outz = std::max(t_zmax, t_zmin);

	// check if t_inx > t_outz
	// OR check if t_inz > t_outx
	// return false, bc missed hit
	if ((t_inx > t_outz) || (t_inz > t_outx)) {
		return false;
	}

	// checks if t_inz > t_inx
	// if true, t_inx = t_inz, bc new lowest t_in
	if (t_inz > t_inx) {
		t_inx = t_inz;
	}

	// checks if t_outz > t_outx
	// if true, t_outx = t_outz, bc new highest t_out
	if (t_outx > t_outz) {
		t_outx = t_outz;
	}

	// hit!
	return true;
}
