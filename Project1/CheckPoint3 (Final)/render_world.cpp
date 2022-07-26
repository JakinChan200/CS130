#include "render_world.h"
#include "flat_shader.h"
#include "object.h"
#include "light.h"
#include "ray.h"

extern bool disable_hierarchy;

Render_World::Render_World()
    :background_shader(0),ambient_intensity(0),enable_shadows(true),
    recursion_depth_limit(3)
{}

Render_World::~Render_World()
{
    delete background_shader;
    for(size_t i=0;i<objects.size();i++) delete objects[i];
    for(size_t i=0;i<lights.size();i++) delete lights[i];
}

// Find and return the Hit structure for the closest intersection.  Be careful
// to ensure that hit.dist>=small_t.
Hit Render_World::Closest_Intersection(const Ray& ray)
{
    //notes 7/27
    Hit o = {nullptr, 0, 0};
    double currentClosestDist = std::numeric_limits<double>::max();

    //For every object in objects, see if ray intersects and if the intersection distance is less than what we currently have
    //If both conditions are true, update o to be the current closest intersection
    for(unsigned long int i = 0; i < objects.size(); i++){
        Hit temp = objects[i]->Intersection(ray, -1);
        if(temp.dist >= small_t && temp.dist < currentClosestDist){
            o = temp;
            currentClosestDist = temp.dist;
        }
    }
    return o;
}

// set up the initial view ray and call
void Render_World::Render_Pixel(const ivec2& pixel_index)
{
    Ray ray;
    //set the endpoint as the camera position
    ray.endpoint = camera.position;
    //set the direction as the normalized path to the pixel
    ray.direction = (camera.World_Position(pixel_index)-ray.endpoint).normalized();
    vec3 color=Cast_Ray(ray,1);
    camera.Set_Pixel(pixel_index,Pixel_Color(color));
}

void Render_World::Render()
{
    if(!disable_hierarchy)
        Initialize_Hierarchy();

    for(int j=0;j<camera.number_pixels[1];j++)
        for(int i=0;i<camera.number_pixels[0];i++)
            Render_Pixel(ivec2(i,j));
}

// cast ray and return the color of the closest intersected surface point,
// or the background color if there is no object intersection
vec3 Render_World::Cast_Ray(const Ray& ray,int recursion_depth)
{
    vec3 color;
    Hit hit = Closest_Intersection(ray);

    //If an intersection existed
    if(hit.object != nullptr){
        //Shade the closest intersection
        vec3 point = (hit.dist * ray.direction) + ray.endpoint;
        color = hit.object->material_shader->Shade_Surface(ray, point, hit.object->Normal(point, hit.part), recursion_depth);
    }else{
        //Shade it based on the background
        color = background_shader->Shade_Surface(ray, ray.direction, ray.direction, recursion_depth);
    }
    return color;
}

void Render_World::Initialize_Hierarchy()
{
    for(auto object : objects){
        for(int i = 0; i < object->number_parts; i++){
            Entry temp;
            temp.obj = object;
            temp.box = object->Bounding_Box(i);
            temp.part = i;
            hierarchy.entries.push_back(temp);
        }
    }

    hierarchy.Reorder_Entries();
    hierarchy.Build_Tree();
}
