//--------------------------
//
// TP done by Arnault VERRET
//--------------------------

#ifndef MESH_H
#define MESH_H

#include <glad/glad.h>
#include <vector>
#include <memory>

#include <glm/glm.hpp>
#include <glm/ext.hpp>

#include <map>
#include <set>

#include <math.h>
#include <iostream>

class Mesh {
public:
  virtual ~Mesh();

  const std::vector<glm::vec3> &vertexPositions() const { return _vertexPositions; }
  std::vector<glm::vec3> &vertexPositions() { return _vertexPositions; }

  const std::vector<glm::vec3> &vertexNormals() const { return _vertexNormals; }
  std::vector<glm::vec3> &vertexNormals() { return _vertexNormals; }

  const std::vector<glm::vec2> &vertexTexCoords() const { return _vertexTexCoords; }
  std::vector<glm::vec2> &vertexTexCoords() { return _vertexTexCoords; }

  const std::vector<glm::uvec3> &triangleIndices() const { return _triangleIndices; }
  std::vector<glm::uvec3> &triangleIndices() { return _triangleIndices; }

  /// Compute the parameters of a sphere which bounds the mesh
  void computeBoundingSphere(glm::vec3 &center, float &radius) const;

  void recomputePerVertexNormals(bool angleBased = false);
  void recomputePerVertexTextureCoordinates( );

  void init();
  void initOldGL();
  void render();
  void clear();

  void addPlan(float square_half_side = 1.0f);

  void subdivideLinear() {
    std::vector<glm::vec3> newVertices = _vertexPositions;
    std::vector<glm::uvec3> newTriangles;

    struct Edge {
      unsigned int a , b;
      Edge( unsigned int c , unsigned int d ) : a( std::min<unsigned int>(c,d) ) , b( std::max<unsigned int>(c,d) ) {}
      bool operator < ( Edge const & o ) const {   return a < o.a  ||  (a == o.a && b < o.b);  }
      bool operator == ( Edge const & o ) const {   return a == o.a  &&  b == o.b;  }
    };
    std::map< Edge , unsigned int > newVertexOnEdge;
    for(unsigned int tIt = 0 ; tIt < _triangleIndices.size() ; ++tIt) {
      unsigned int a = _triangleIndices[tIt][0];
      unsigned int b = _triangleIndices[tIt][1];
      unsigned int c = _triangleIndices[tIt][2];


      Edge Eab(a,b);
      unsigned int oddVertexOnEdgeEab = 0;
      if( newVertexOnEdge.find( Eab ) == newVertexOnEdge.end() ) {
        newVertices.push_back( (_vertexPositions[ a ] + _vertexPositions[ b ]) / 2.f );
        oddVertexOnEdgeEab = newVertices.size() - 1;
        newVertexOnEdge[Eab] = oddVertexOnEdgeEab;
      }
      else { oddVertexOnEdgeEab = newVertexOnEdge[Eab]; }


      Edge Ebc(b,c);
      unsigned int oddVertexOnEdgeEbc = 0;
      if( newVertexOnEdge.find( Ebc ) == newVertexOnEdge.end() ) {
        newVertices.push_back( (_vertexPositions[ b ] + _vertexPositions[ c ]) / 2.f );
        oddVertexOnEdgeEbc = newVertices.size() - 1;
        newVertexOnEdge[Ebc] = oddVertexOnEdgeEbc;
      }
      else { oddVertexOnEdgeEbc = newVertexOnEdge[Ebc]; }

      Edge Eca(c,a);
      unsigned int oddVertexOnEdgeEca = 0;
      if( newVertexOnEdge.find( Eca ) == newVertexOnEdge.end() ) {
        newVertices.push_back( (_vertexPositions[ c ] + _vertexPositions[ a ]) / 2.f );
        oddVertexOnEdgeEca = newVertices.size() - 1;
        newVertexOnEdge[Eca] = oddVertexOnEdgeEca;
      }
      else { oddVertexOnEdgeEca = newVertexOnEdge[Eca]; }

      // set new triangles :
      newTriangles.push_back( glm::uvec3( a , oddVertexOnEdgeEab , oddVertexOnEdgeEca ) );
      newTriangles.push_back( glm::uvec3( oddVertexOnEdgeEab , b , oddVertexOnEdgeEbc ) );
      newTriangles.push_back( glm::uvec3( oddVertexOnEdgeEca , oddVertexOnEdgeEbc , c ) );
      newTriangles.push_back( glm::uvec3( oddVertexOnEdgeEab , oddVertexOnEdgeEbc , oddVertexOnEdgeEca ) );
    }

    // after that:
    _triangleIndices = newTriangles;
    _vertexPositions = newVertices;
    recomputePerVertexNormals( );
    recomputePerVertexTextureCoordinates( );
  }

  /* Loop subdivision method */
  void subdivideLoop(){
    std::vector<glm::vec3> newVertices = _vertexPositions;
    std::vector<glm::uvec3> newTriangles;

    struct Edge {
      unsigned int a , b;
      Edge( unsigned int c , unsigned int d ) : a( std::min<unsigned int>(c,d) ) , b( std::max<unsigned int>(c,d) ) {}
      bool operator < ( Edge const & o ) const {   return a < o.a  ||  (a == o.a && b < o.b);  }
      bool operator == ( Edge const & o ) const {   return a == o.a  &&  b == o.b;  }
    };

    // This variable stores all the edges of each even vertices
    // For each edges associated, we store how many times we saw it in the mesh
    // Thus, we can easily check if the vertex is ordinary or not
    std::vector<std::map<Edge, unsigned int>> evenEdgeNeighbors(_vertexPositions.size());

    // We loop to fill evenEdgeNeighbors
    for(unsigned int tIt = 0 ; tIt < _triangleIndices.size() ; ++tIt) {
      unsigned int a = _triangleIndices[tIt][0];
      unsigned int b = _triangleIndices[tIt][1];
      unsigned int c = _triangleIndices[tIt][2];
      Edge Eab(a, b);
      Edge Ebc(b, c);
      Edge Eac(a, c);

      // a
      evenEdgeNeighbors[a][Eab] = (evenEdgeNeighbors[a].find(Eab) == evenEdgeNeighbors[a].end()) ? 1 : 2;
      evenEdgeNeighbors[a][Eac] = (evenEdgeNeighbors[a].find(Eac) == evenEdgeNeighbors[a].end()) ? 1 : 2;
      // b
      evenEdgeNeighbors[b][Eab] = (evenEdgeNeighbors[b].find(Eab) == evenEdgeNeighbors[b].end()) ? 1 : 2;
      evenEdgeNeighbors[b][Ebc] = (evenEdgeNeighbors[b].find(Ebc) == evenEdgeNeighbors[b].end()) ? 1 : 2;
      // c
      evenEdgeNeighbors[c][Eac] = (evenEdgeNeighbors[c].find(Eac) == evenEdgeNeighbors[c].end()) ? 1 : 2;
      evenEdgeNeighbors[c][Ebc] = (evenEdgeNeighbors[c].find(Ebc) == evenEdgeNeighbors[c].end()) ? 1 : 2;
    }

    // we move the even vertices
    for (unsigned int i = 0; i < _vertexPositions.size(); i++) {
      bool is_ordinary = true;
      int n = evenEdgeNeighbors[i].size(); // number of neighbors
      float a_n = (40.f - pow(3.f + 2.f*cos(2.f*M_PI/n), 2)) / 64.f;
      newVertices[i] *= (1 - a_n);


      // We'll consider the extraordinary case only if we face an edge appearing once.
      // Thus, we can just assume it is ordinary and then scrap everything we've done
      // if we face an edge with a value of 1.
      // as a_n is already computed, it doesn't cost so much as we'd have to
      // loop on all edges anyway.
      for (auto it = evenEdgeNeighbors[i].begin(); it != evenEdgeNeighbors[i].end(); ++it) {
        unsigned int neighbor_vertex = (it->first.a == i) ? it->first.b : it->first.a;
        if (is_ordinary){ // everything has been ok at this point
          if (it->second == 2){
            newVertices[i] += _vertexPositions[neighbor_vertex]*a_n/(float)n;
          }
          else{ // we face the first edge appearing once
            is_ordinary = false;
            newVertices[i] = _vertexPositions[i]*3.f/4.f + _vertexPositions[neighbor_vertex]/8.f;
          }
        }
        else{ // we got a solo edge, we assume there can be only 2 edges appearing once.
          if(it->second != 2)
            newVertices[i] += _vertexPositions[neighbor_vertex]/8.f;
        }
      }
    }


    // calculate odd vertices
    // to know if an odd vertex is ordinary, we keep track of the vertex opposite to the edge of the odd vertex
    // Then, if we pass again through the edge, we can take in count both opposite vertices
    std::map< Edge , unsigned int > newVertexOnEdge;
    std::map<unsigned int, unsigned int> oddValence;
    for(unsigned int tIt = 0 ; tIt < _triangleIndices.size() ; ++tIt) {
      unsigned int a = _triangleIndices[tIt][0];
      unsigned int b = _triangleIndices[tIt][1];
      unsigned int c = _triangleIndices[tIt][2];

      // We consider each edges
      Edge Eab(a,b);
      unsigned int oddVertexOnEdgeEab = 0;
      if( newVertexOnEdge.find( Eab ) == newVertexOnEdge.end() ) {
        newVertices.push_back( (_vertexPositions[ a ] + _vertexPositions[ b ]) / 2.f );
        oddVertexOnEdgeEab = newVertices.size() - 1;
        newVertexOnEdge[Eab] = oddVertexOnEdgeEab;
        oddValence[oddVertexOnEdgeEab] = c; // keep track of the 3rd vertex in case of ordinary odd vertex
      }
      else {
        oddVertexOnEdgeEab = newVertexOnEdge[Eab];
        newVertices[oddVertexOnEdgeEab] = newVertices[oddVertexOnEdgeEab]*3.f/4.f + // 1/2 * 3/4 = 3/8
                                          _vertexPositions[oddValence[oddVertexOnEdgeEab]]/8.f +
                                          _vertexPositions[c]/8.f; //change value for ordinary value
      }


      Edge Ebc(b,c);
      unsigned int oddVertexOnEdgeEbc = 0;
      if( newVertexOnEdge.find( Ebc ) == newVertexOnEdge.end() ) {
        newVertices.push_back( (_vertexPositions[ b ] + _vertexPositions[ c ]) / 2.f );
        oddVertexOnEdgeEbc = newVertices.size() - 1;
        newVertexOnEdge[Ebc] = oddVertexOnEdgeEbc;
        oddValence[oddVertexOnEdgeEbc] = a; // keep track of the 3rd vertex in case of ordinary odd vertex

      }
      else { oddVertexOnEdgeEbc = newVertexOnEdge[Ebc];
      newVertices[oddVertexOnEdgeEbc] = newVertices[oddVertexOnEdgeEbc]*3.f/4.f +
                                        _vertexPositions[oddValence[oddVertexOnEdgeEbc]]/8.f +
                                        _vertexPositions[a]/8.f; //change value for ordinary value
      }

      Edge Eca(c,a);
      unsigned int oddVertexOnEdgeEca = 0;
      if( newVertexOnEdge.find( Eca ) == newVertexOnEdge.end() ) {
        newVertices.push_back( (_vertexPositions[ c ] + _vertexPositions[ a ]) / 2.f );
        oddVertexOnEdgeEca = newVertices.size() - 1;
        newVertexOnEdge[Eca] = oddVertexOnEdgeEca;
        oddValence[oddVertexOnEdgeEca] = b; // keep track of the 3rd vertex in case of ordinary odd vertex
      }
      else { oddVertexOnEdgeEca = newVertexOnEdge[Eca];
      newVertices[oddVertexOnEdgeEca] = newVertices[oddVertexOnEdgeEca]*3.f/4.f +
                                        _vertexPositions[oddValence[oddVertexOnEdgeEca]]/8.f +
                                        _vertexPositions[b]/8.f; //change value for ordinary value
      }

      // set new triangles :
      newTriangles.push_back( glm::uvec3( a , oddVertexOnEdgeEab , oddVertexOnEdgeEca ) );
      newTriangles.push_back( glm::uvec3( oddVertexOnEdgeEab , b , oddVertexOnEdgeEbc ) );
      newTriangles.push_back( glm::uvec3( oddVertexOnEdgeEca , oddVertexOnEdgeEbc , c ) );
      newTriangles.push_back( glm::uvec3( oddVertexOnEdgeEab , oddVertexOnEdgeEbc , oddVertexOnEdgeEca ) );
    }

    // after that:
    _triangleIndices = newTriangles;
    _vertexPositions = newVertices;
    recomputePerVertexNormals( );
    recomputePerVertexTextureCoordinates( );
  }

  /* Catmull Clark subdivision method */
  void subdivideCatmullClark(){
    std::vector<glm::vec3> newVertices = _vertexPositions;
    std::vector<glm::uvec3> newTriangles;

    struct Edge {
      unsigned int a , b;
      Edge( unsigned int c , unsigned int d ) : a( std::min<unsigned int>(c,d) ) , b( std::max<unsigned int>(c,d) ) {}
      bool operator < ( Edge const & o ) const {   return a < o.a  ||  (a == o.a && b < o.b);  }
      bool operator == ( Edge const & o ) const {   return a == o.a  &&  b == o.b;  }
    };

    // face points
    // we keep track of vertices linked to the face point
    std::vector<std::vector<unsigned int>> evenFaceNeighbors(_vertexPositions.size());
    std::map< Edge , unsigned int > newVertexOnEdge;
    std::map<unsigned int, unsigned int> oddValence;
    for(unsigned int tIt = 0 ; tIt < _triangleIndices.size() ; ++tIt) {
      unsigned int a = _triangleIndices[tIt][0];
      unsigned int b = _triangleIndices[tIt][1];
      unsigned int c = _triangleIndices[tIt][2];

      newVertices.push_back((_vertexPositions[a]+_vertexPositions[b]+_vertexPositions[c])/3.f);
      unsigned int face_id = newVertices.size() - 1;
      evenFaceNeighbors[a].push_back(face_id);
      evenFaceNeighbors[b].push_back(face_id);
      evenFaceNeighbors[c].push_back(face_id);


      // We consider each edges
      Edge Eab(a,b);
      unsigned int oddVertexOnEdgeEab = 0;
      if( newVertexOnEdge.find( Eab ) == newVertexOnEdge.end() ) {
        newVertices.push_back( (_vertexPositions[ a ] + _vertexPositions[ b ]) / 2.f );
        oddVertexOnEdgeEab = newVertices.size() - 1;
        newVertexOnEdge[Eab] = oddVertexOnEdgeEab;
        oddValence[oddVertexOnEdgeEab] = face_id; // keep track of the face vertex in case of ordinary odd vertex
      }
      else {
        oddVertexOnEdgeEab = newVertexOnEdge[Eab];
        newVertices[oddVertexOnEdgeEab] = newVertices[oddVertexOnEdgeEab]/2.f +
                                          newVertices[oddValence[oddVertexOnEdgeEab]]/4.f +
                                          newVertices[face_id]/4.f; //change value for ordinary value
      }


      Edge Ebc(b,c);
      unsigned int oddVertexOnEdgeEbc = 0;
      if( newVertexOnEdge.find( Ebc ) == newVertexOnEdge.end() ) {
        newVertices.push_back( (_vertexPositions[ b ] + _vertexPositions[ c ]) / 2.f );
        oddVertexOnEdgeEbc = newVertices.size() - 1;
        newVertexOnEdge[Ebc] = oddVertexOnEdgeEbc;
        oddValence[oddVertexOnEdgeEbc] = face_id; // keep track of the face vertex in case of ordinary odd vertex

      }
      else { oddVertexOnEdgeEbc = newVertexOnEdge[Ebc];
      newVertices[oddVertexOnEdgeEbc] = newVertices[oddVertexOnEdgeEbc]/2.f +
                                        newVertices[oddValence[oddVertexOnEdgeEbc]]/4.f +
                                        newVertices[face_id]/4.f; //change value for ordinary value
      }

      Edge Eca(c,a);
      unsigned int oddVertexOnEdgeEca = 0;
      if( newVertexOnEdge.find( Eca ) == newVertexOnEdge.end() ) {
        newVertices.push_back( (_vertexPositions[ c ] + _vertexPositions[ a ]) / 2.f );
        oddVertexOnEdgeEca = newVertices.size() - 1;
        newVertexOnEdge[Eca] = oddVertexOnEdgeEca;
        oddValence[oddVertexOnEdgeEca] = face_id; // keep track of the face vertex in case of ordinary odd vertex
      }
      else { oddVertexOnEdgeEca = newVertexOnEdge[Eca];
      newVertices[oddVertexOnEdgeEca] = newVertices[oddVertexOnEdgeEca]/2.f +
                                        newVertices[oddValence[oddVertexOnEdgeEca]]/4.f +
                                        newVertices[face_id]/4.f; //change value for ordinary value
      }


      // set new triangles :
      newTriangles.push_back( glm::uvec3( a , oddVertexOnEdgeEab , face_id ) );
      newTriangles.push_back( glm::uvec3( oddVertexOnEdgeEab , b , face_id ) );
      newTriangles.push_back( glm::uvec3( face_id , b , oddVertexOnEdgeEbc ) );
      newTriangles.push_back( glm::uvec3( a , face_id , oddVertexOnEdgeEca ) );
      newTriangles.push_back( glm::uvec3( c , face_id , oddVertexOnEdgeEbc ) );
      newTriangles.push_back( glm::uvec3( face_id , c , oddVertexOnEdgeEca ) );
    }

    //move even points
    // This variable stores all the edges of each even vertices
    // For each edges associated, we store how many times we saw it in the mesh
    // Thus, we can easily check if the vertex is ordinary or not
    std::vector<std::map<Edge, unsigned int>> evenEdgeNeighbors(_vertexPositions.size());

    // We loop to fill evenEdgeNeighbors
    for(unsigned int tIt = 0 ; tIt < _triangleIndices.size() ; ++tIt) {
      unsigned int a = _triangleIndices[tIt][0];
      unsigned int b = _triangleIndices[tIt][1];
      unsigned int c = _triangleIndices[tIt][2];
      Edge Eab(a, b);
      Edge Ebc(b, c);
      Edge Eac(a, c);

      // a
      evenEdgeNeighbors[a][Eab] = (evenEdgeNeighbors[a].find(Eab) == evenEdgeNeighbors[a].end()) ? 1 : 2;
      evenEdgeNeighbors[a][Eac] = (evenEdgeNeighbors[a].find(Eac) == evenEdgeNeighbors[a].end()) ? 1 : 2;
      // b
      evenEdgeNeighbors[b][Eab] = (evenEdgeNeighbors[b].find(Eab) == evenEdgeNeighbors[b].end()) ? 1 : 2;
      evenEdgeNeighbors[b][Ebc] = (evenEdgeNeighbors[b].find(Ebc) == evenEdgeNeighbors[b].end()) ? 1 : 2;
      // c
      evenEdgeNeighbors[c][Eac] = (evenEdgeNeighbors[c].find(Eac) == evenEdgeNeighbors[c].end()) ? 1 : 2;
      evenEdgeNeighbors[c][Ebc] = (evenEdgeNeighbors[c].find(Ebc) == evenEdgeNeighbors[c].end()) ? 1 : 2;
    }

    // we move the even vertices
    for (unsigned int i = 0; i < _vertexPositions.size(); i++) {
      bool is_ordinary = true;
      int n = evenEdgeNeighbors[i].size(); // number of neighbors
      newVertices[i] = (float)(n-2)*_vertexPositions[i]/(float)n;

      // We'll consider the extraordinary case only if we face an edge appearing once.
      // Thus, we can just assume it is ordinary and then scrap everything we've done
      // if we face an edge with a value of 1.
      for (auto it = evenEdgeNeighbors[i].begin(); it != evenEdgeNeighbors[i].end(); ++it) {
        unsigned int neighbor_vertex = (it->first.a == i) ? it->first.b : it->first.a;
        if (is_ordinary){ // everything has been ok at this point
          if (it->second == 2){
            newVertices[i] += _vertexPositions[neighbor_vertex]/(float)evenEdgeNeighbors[i].size()/(float)n;
          }
          else{ // we face the first edge appearing once
            is_ordinary = false;
            newVertices[i] = _vertexPositions[i]*3.f/4.f + _vertexPositions[neighbor_vertex]/8.f;
          }
        }
        else{ // we got a solo edge, we assume there can be only 2 edges appearing once.
          if(it->second != 2)
            newVertices[i] += _vertexPositions[neighbor_vertex]/8.f;
        }
      }

      if(is_ordinary){
        for (unsigned int j = 0; j < evenFaceNeighbors[i].size(); ++j) {
          newVertices[i] += newVertices[evenFaceNeighbors[i][j]]/(float)evenFaceNeighbors[i].size()/(float)n;
        }
      }
    }
    // after that:
    _triangleIndices = newTriangles;
    _vertexPositions = newVertices;
    recomputePerVertexNormals( );
    recomputePerVertexTextureCoordinates( );
  }

  void subdivide() {
    //subdivideLinear();
    subdivideLoop();
    //subdivideCatmullClark();
  }

private:
  std::vector<glm::vec3> _vertexPositions;
  std::vector<glm::vec3> _vertexNormals;
  std::vector<glm::vec2> _vertexTexCoords;
  std::vector<glm::uvec3> _triangleIndices;

  GLuint _vao = 0;
  GLuint _posVbo = 0;
  GLuint _normalVbo = 0;
  GLuint _texCoordVbo = 0;
  GLuint _ibo = 0;
};

// utility: loader
void loadOFF(const std::string &filename, std::shared_ptr<Mesh> meshPtr);

#endif  // MESH_H
