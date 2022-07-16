#include "light.h"
#include "object.h"
#include "phong_shader.h"
#include "ray.h"
#include "render_world.h"

vec3 Phong_Shader::
Shade_Surface(const Ray& ray,const vec3& intersection_point,
    const vec3& normal,int recursion_depth) const
{
    vec3 color;
    
    //color = RL + RLmax(dot(n, l), 0) + RLmax(dot(v, r), 0)^a
    //L = light intensity
    //R = reflectivity

    color = world.ambient_intensity * world.ambient_color * color_ambient;

    for(long unsigned int i = 0 ; i < world.lights.size(); i++){
        vec3 l = (world.lights[i]->position - intersection_point);

        double diffuseScalar = std::max(dot(normal, l.normalized()), 0.0);
        color += world.lights[i]->Emitted_Light(l) * color_diffuse * diffuseScalar;
        
    }

    //color = color_ambient + color_diffuse + color_specular;
    //TODO; //determine the color
    //Ambient + diffuse + Specular = phong
    //Ambient: I = LaRa
    return color;
}
