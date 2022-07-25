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

        if(world.enable_shadows){
            Ray shadowRay(intersection_point, l);
            Hit shadowIntersect = world.Closest_Intersection(shadowRay);
            if(shadowIntersect.object != nullptr && shadowIntersect.dist < l.magnitude()){
                continue;
            }
        }

        double diffuseScalar = std::max(dot(normal, l.normalized()), 0.0);
        color += world.lights[i]->Emitted_Light(l) * color_diffuse * diffuseScalar;

        vec3 r = -(l - (2.0*dot(l, normal) * normal)).normalized();
        vec3 v = (world.camera.position - intersection_point).normalized();
        double specularScalar = pow(std::max(dot(v, r), 0.0), specular_power);
        color += world.lights[i]->Emitted_Light(l) * specularScalar * color_specular;
    }

    //color_specular
    //specular_power
    //r = l - 2dot(l, n)/dot(n, n) * n      can cross out denominator if normal is unit vector https://www.youtube.com/watch?v=naaeH1qbjdQ&t=650s

    //color = color_ambient + color_diffuse + color_specular;
    //TODO; //determine the color
    //Ambient + diffuse + Specular = phong
    //Ambient: I = LaRa
    return color;
}
