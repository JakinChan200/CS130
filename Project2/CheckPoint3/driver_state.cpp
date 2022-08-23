#include "driver_state.h"
#include <cstring>
#include <iostream>
#include <vector>
#include <limits>

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
    state.image_depth= new float[dimensions];

    for(int i = 0; i < dimensions; i++){
        state.image_color[i] = make_pixel(0, 0, 0);
        state.image_depth[i] = std::numeric_limits<float>::max();
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
    int vert = 0;

    switch(type){
        case render_type::triangle:{
            for(int i = 0; i < state.num_vertices / 3; i++){
                for(int j = 0; j < 3; j++){
                    custData[j].data = &state.vertex_data[(i * 3 + j)* state.floats_per_vertex];
                    temp[j].data = custData[j].data;
                    state.vertex_shader(custData[j], temp[j], state.uniform_data);
                }
                clip_triangle(state, temp[0], temp[1], temp[2], 0);            
            }
            break;
        }
        case render_type::indexed:{ //index_data holds the index for each vertex info from vertex_data
            for(int i = 0; i < 3*state.num_triangles; i++){
                vert = i % 3;
                custData[vert].data = &state.vertex_data[state.index_data[i]* state.floats_per_vertex];
                temp[vert].data = custData[vert].data;
                state.vertex_shader(custData[vert], temp[vert], state.uniform_data);
                if(i != 0 && vert == 2){ //If all three vertices are initialized full (technically doenst check for it though)
                    clip_triangle(state, temp[0], temp[1], temp[2], 0);
                }          
            }
            break;
        }
        case render_type::fan:{ //Keep the first vertex as the common first one for all triangles
            for(int i = 0; i < state.num_vertices; i++){
                for(int j = 0; j < 3; j++){
                    custData[j].data = j == 0 ? state.vertex_data : &state.vertex_data[((i + j) * state.floats_per_vertex)];
                    temp[j].data = custData[j].data;
                    state.vertex_shader(custData[j], temp[j], state.uniform_data);
                }
                clip_triangle(state, temp[0], temp[1], temp[2], 0);
            }
            break;
        }
        case render_type::strip:{   //Instead of 3N, its N + 2 where N is the number of vertices (ABCDE -> ABC, BCD, CDE)
            for(int i = 0; i < state.num_vertices-2; i++){
                for(int j = 0; j < 3; j++){
                    custData[j].data = &state.vertex_data[(i+j)* state.floats_per_vertex];
                    temp[j].data = custData[j].data;
                    state.vertex_shader(custData[j], temp[j], state.uniform_data);
                }
                clip_triangle(state, temp[0], temp[1], temp[2], 0);
            }
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
    const data_geometry& v1, const data_geometry& v2,int face){
    
    //https://www.ijirt.org/master/publishedpaper/IJIRT100119_PAPER.pdf
    if(face==6){
        rasterize_triangle(state, v0, v1, v2);
        return;
    }

    vec4 A = v0.gl_Position;
    vec4 B = v1.gl_Position;
    vec4 C = v2.gl_Position;

    //Innocent till proven guilty
    bool isAIn = true;
    bool isBIn = true;
    bool isCIn = true;
    float sign = 1;
    int axis = 0;

    switch(face){
        case 0: // x <= w       right side
            isAIn = A[0] <= A[3] ? true : false;
            isBIn = B[0] <= B[3] ? true : false;
            isCIn = C[0] <= C[3] ? true : false;
            axis = 0;
            break;
        case 1: // x >= -w      left side
            isAIn = A[0] >= -A[3] ? true : false;
            isBIn = B[0] >= -B[3] ? true : false;
            isCIn = C[0] >= -C[3] ? true : false;
            sign = -1;
            axis = 0;
            break;
        case 2: // y <= w       top side
            isAIn = A[1] <= A[3] ? true : false;
            isBIn = B[1] <= B[3] ? true : false;
            isCIn = C[1] <= C[3] ? true : false;
            axis = 1;
            break;
        case 3: // y >= -w      bot side
            isAIn = A[1] >= -A[3] ? true : false;
            isBIn = B[1] >= -B[3] ? true : false;
            isCIn = C[1] >= -C[3] ? true : false;
            sign = -1;
            axis = 1;
            break;
        case 4: // z <= w       far side
            isAIn = A[2] <= A[3] ? true : false;
            isBIn = B[2] <= B[3] ? true : false;
            isCIn = C[2] <= C[3] ? true : false;
            axis = 2;
            break;
        case 5: //z >= -w       near side
            isAIn = A[2] >= -A[3] ? true : false;
            isBIn = B[2] >= -B[3] ? true : false;
            isCIn = C[2] >= -C[3] ? true : false;
            sign = -1;
            axis = 2;  
            break;
        default:
            break;
    }

    /*
    I I I
    I I 0
    I 0 I
    I 0 0
    0 I I
    0 I 0
    0 0 I
    0 0 0
    */
    data_geometry vert1;
    data_geometry vert2;

    if(isAIn && isBIn && isCIn){
        clip_triangle(state, v0, v1, v2, face+1);
    }else if(isAIn && isBIn && !isCIn){
        vert1 = createTriangle(state, v0, v2, axis, sign);
        vert2 = createTriangle(state, v1, v2, axis, sign);
        clip_triangle(state, v0, v1, vert1, face+1);
        clip_triangle(state, v1, vert1, vert2, face+1);
    }else if(isAIn && !isBIn && isCIn){
        vert1 = createTriangle(state, v0, v1, axis, sign);
        vert2 = createTriangle(state, v2, v1, axis, sign);
        clip_triangle(state, v0, vert1, v2, face+1);
        clip_triangle(state, vert1, vert2, v2, face+1);
    }else if(isAIn && !isBIn && !isCIn){
        vert1 = createTriangle(state, v0, v1, axis, sign);
        vert2 = createTriangle(state, v0, v2, axis, sign);
        clip_triangle(state, v0, vert1, vert2, face+1);
    }else if(!isAIn && isBIn && isCIn){
        vert1 = createTriangle(state, v2, v0, axis, sign);
        vert2 = createTriangle(state, v1, v0, axis, sign);
        clip_triangle(state, vert2, v1, v2, face+1);
        clip_triangle(state, vert2, v2, vert1, face+1);
    }else if(!isAIn && isBIn && !isCIn){
        vert1 = createTriangle(state, v1, v0, axis, sign);
        vert2 = createTriangle(state, v1, v2, axis, sign);
        clip_triangle(state, vert1, v1, vert2, face+1);
    }else if(!isAIn && !isBIn && isCIn){
        vert1 = createTriangle(state, v2, v0, axis, sign);
        vert2 = createTriangle(state, v2, v1, axis, sign);
        clip_triangle(state, vert1, vert2, v2, face+1);
    }else if(!isAIn && !isBIn && !isCIn){
        return;
    }else{
        return;
    }

    //std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    //clip_triangle(state,v0,v1,v2,face+1);
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
    float alpha, beta, gamma;

    auto *data = new float[MAX_FLOATS_PER_VERTEX];
    data_fragment pixelData{data};
    data_output color;
    float zBuffer;
    float denom;

    for(int i = minY; i <= maxY; i++){
        for(int j = minX; j <= maxX; j++){
            alpha = 0.5f*(((v12[0] * v22[1]) - (v22[0] * v12[1])) + ((v22[0] * i) - (j * v22[1])) + ((j*v12[1]) - (v12[0]*i)))/triangleArea;
            beta = 0.5f*(((v22[0] * v02[1]) - (v02[0] * v22[1])) + ((v02[0] * i) - (j * v02[1])) + ((j*v22[1]) - (v22[0]*i)))/triangleArea;
            gamma = 0.5f*(((v02[0] * v12[1]) - (v12[0] * v02[1])) + ((v12[0] * i) - (j * v12[1])) + ((j*v02[1]) - (v02[0]*i)))/triangleArea;
            
            if(alpha >= 0 && beta >= 0 && gamma >= 0){
                zBuffer = alpha*v0.gl_Position[2]/v0.gl_Position[3] + beta*v1.gl_Position[2]/v1.gl_Position[3] + gamma*v2.gl_Position[2]/v2.gl_Position[3]; //Only take the closest color
                if(!(zBuffer < state.image_depth[j + (i*state.image_width)])){continue;}

                for(int k = 0; k < state.floats_per_vertex; k++){
                    switch(state.interp_rules[k]){
                        case interp_type::flat:
                            pixelData.data[k] = v0.data[k]; //Use the data of the first vertex
                            break;
                        case interp_type::smooth:
                            denom = alpha/v0.gl_Position[3] + beta/v1.gl_Position[3] + gamma/v2.gl_Position[3];
                            
                            pixelData.data[k] = ((alpha / v0.gl_Position[3] / denom) * v0.data[k]) 
                                              + ((beta  / v1.gl_Position[3] / denom) * v1.data[k]) 
                                              + ((gamma / v2.gl_Position[3] / denom) * v2.data[k]);
                            break;
                        case interp_type::noperspective:
                            pixelData.data[k] = alpha*v0.data[k] + beta*v1.data[k] + gamma*v2.data[k]; //Use the data of all three, with no perspective changes
                            break;
                        default:
                            break;
                    }
                }
                state.fragment_shader(pixelData, color, state.uniform_data);

                state.image_depth[j + (i*state.image_width)] = zBuffer;
                state.image_color[j + (i*state.image_width)] = make_pixel(color.output_color[0]*255, color.output_color[1]*255, color.output_color[2]*255);
            }
        }
    }
    //std::cout<<"TODO: implement rasterization"<<std::endl;
}

data_geometry createTriangle(driver_state& state, const data_geometry& v0, const data_geometry& v1, int axis, float sign){
    data_geometry temp;
    auto *data = new float[MAX_FLOATS_PER_VERTEX];

    //https://youtu.be/VMD7fsCYO9o?t=854
    float alpha = ((sign * v0.gl_Position[3]) - v0.gl_Position[axis]) / ((sign * v0.gl_Position[3]) - v0.gl_Position[axis] - (sign * v1.gl_Position[3]) + v1.gl_Position[axis]);

    //https://www.cs.ucr.edu/~craigs/courses/2022-summer-cs-130/lectures/barycentric-coordinates.pdf
    temp.gl_Position = alpha * v1.gl_Position + (1 - alpha) * v0.gl_Position;
    float noPersAlpha = (alpha * v1.gl_Position[3]) / (alpha * v1.gl_Position[3] + (1 - alpha) * v0.gl_Position[3]);

    for (int i = 0; i < state.floats_per_vertex; i++) {
		switch(state.interp_rules[i]) {
            case(interp_type::flat):
                data[i] = v0.data[i];
                break;
            case(interp_type::smooth):
                data[i] = alpha * v1.data[i] + (1 - alpha) * v0.data[i];
                break;
            case(interp_type::noperspective):
                data[i] = noPersAlpha * v1.data[i] + (1 - noPersAlpha) * v0.data[i];
                break;
            default:
                break;
        }
	}
    temp.data = data;
    return temp;
}