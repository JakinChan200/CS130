#include "driver_state.h"
#include <cstring>
#include <iostream>
#include <vector>

using namespace std;
driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width=width;
    state.image_height=height;

    int dimensions = width * height;
    state.image_color= new pixel[dimensions];
    state.image_depth= 0;//new float[dimensions];

    for(int i = 0; i < dimensions; i++){
        state.image_color[i] = make_pixel(0, 0, 0);
    }

    
    //std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
    data_geometry temp[3];
    data_vertex custData[3];

    switch(type){
        case render_type::triangle:{
            for (int i = 0; i < state.num_vertices / 3; i++) {
                for (int j = 0; j < 3; j++) {
                    custData[j].data = &state.vertex_data[(i * 3 + j)* state.floats_per_vertex];
                    temp[j].data = custData[j].data;
                    state.vertex_shader(custData[j], temp[j], state.uniform_data);
                }
                rasterize_triangle(state, temp[0], temp[1], temp[2]);
            }
            break;
        }
        case render_type::indexed:{

            break;
        }
        case render_type::fan:{

            break;
        }
        case render_type::strip:{

            break;
        }
        default:
            break;
    }
    //std::cout<<"TODO: implement rendering."<<std::endl;
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2,int face)
{
    if(face==6)
    {
        rasterize_triangle(state, v0, v1, v2);
        return;
    }
    std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state,v0,v1,v2,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2)
{
    float v02[2];
    float v12[2];
    float v22[2];

    float halfWidth = state.image_width/2.0f;
    float halfHeight = state.image_height/2.0f;

    v02[0] = halfWidth * (v0.gl_Position[0]/v0.gl_Position[3]) + halfWidth -0.5f;
    v02[1] = halfHeight * (v0.gl_Position[1]/v0.gl_Position[3]) + halfHeight -0.5f;

    v12[0] = halfWidth * (v1.gl_Position[0]/v1.gl_Position[3]) + halfWidth -0.5f;
    v12[1] = halfHeight * (v1.gl_Position[1]/v1.gl_Position[3]) + halfHeight -0.5f;

    v22[0] = halfWidth * (v2.gl_Position[0]/v2.gl_Position[3]) + halfWidth -0.5f;
    v22[1] = halfHeight * (v2.gl_Position[1]/v2.gl_Position[3]) + halfHeight -0.5f;

    float minX = std::min(std::min(v02[0],v12[0]),v22[0]);
    float maxX = std::max(std::max(v02[0],v12[0]),v22[0]);
    float minY = std::min(std::min(v02[1],v12[1]),v22[1]);
    float maxY = std::max(std::max(v02[1],v12[1]),v22[1]);

    float triangleArea = 0.5f*(((v12[0] * v22[1]) - (v22[0] * v12[1])) + ((v22[0] * v02[1]) - (v02[0] * v22[1])) + ((v02[0]*v12[1]) - (v12[0]*v02[1])));
    float alpha;
    float beta;
    float gamma;

    for(int i = minY; i <= maxY; i++){
        for(int j = minX; j <= maxX; j++){
            alpha = 0.5f*(((v12[0] * v22[1]) - (v22[0] * v12[1])) + ((v22[0] * i) - (j * v22[1])) + ((j*v12[1]) - (v12[0]*i)))/triangleArea;
            beta = 0.5f*(((v22[0] * v02[1]) - (v02[0] * v22[1])) + ((v02[0] * i) - (j * v02[1])) + ((j*v22[1]) - (v22[0]*i)))/triangleArea;
            gamma = 0.5f*(((v02[0] * v12[1]) - (v12[0] * v02[1])) + ((v12[0] * i) - (j * v12[1])) + ((j*v02[1]) - (v02[0]*i)))/triangleArea;
            
            if(alpha >= 0 && beta >= 0 && gamma >= 0){
                state.image_color[j + (i*state.image_width)] = make_pixel(255, 255, 255);
            }
        }
    }
    //std::cout<<"TODO: implement rasterization"<<std::endl;
}

