#ifndef BUCKET_HPP_
#define BUCKET_HPP_

/*
 * @brief Bucket class is a container of entities representing the different components of a mesh
 * The different entites supported are vertices, edges, faces, and elements
 * Each bucket stores information about homogenous entities (i.e, a given bucket can reference both vertices and elements for instance)
 * A bucket stores its capacity which is the maximum number of entites that can be references by the bucket, and is set to a constant defined in the constants.h file
 * The bucket capacity can be changed from the constructor of the bucket
 * A bucket stores its size which is the actual number of entities references by this bucket
 * A bucket holds a pointer to the first entity it references and from which the other entities are references
 * This pointer points to the entity stored in the entity vector that lives in the entity container
 * IT IS IMPORTANT TO UNDERSTAND THAT THE BUCKET CLASS IS CREATED TO GROUP OBJECTS THAT ARE STORED CONTIGUOUSLY IN THE MEMORY OF THE COMPUER
 * DO NOT USE IT ON A LIST OF OBJECTS THAT ARE SPREAD IN THE MEMORY OF THE COMPUTER AS IT WILL RESULT IN UNDETERMINED BEHAVIOR
 */

#include<stdexcept>

#include "Utility/Exceptions.hpp"
#include "Utility/Constants.hpp"
#include "Types/Types.hpp"

namespace NSMesh { namespace NSEntity {

class Entity;

template<typename T>
class Bucket {

protected:
  size_t _size; /**<Number of entities referenced by this bucket*/
  size_t _capacity; /**<Maximum number of entities that can be references by a single bucket*/
  const T* _beginData; /**<Pointer to the first entity referenced by this bucket*/
  std::vector<T> _values = {}; /**<Vector of pointers to the values that will be referenced by this bucket. Note this might not necessarily be needed*/
  
public:
  /*
   * Constructor
   @param contiguous - Parameter specifying if the data to be referenced by the Bucket is contiguous in memory or not
   @param capacity - Maximum number of entities to be referenced by this bucket
  */
  explicit Bucket( size_t p_capacity = BUCKETSIZE );

  /*
   * Function that returns the number of elements referenced by a given bucket
   @output Number of elements referenced by the bucket
  */
  size_t GetSize() const;

  /*
   * Function that returns the maximum number of elements that can be references by a given bucket
   @output Maximum number of elements that can be referenced by the bucket
  */
  size_t GetCapacity() const;

  /*
   * Operator to access the entity of index i within the bucket
   @param i - index of the entity to be accessed
   @output The elemenet at index i in the bucket by reference
  */
  const T& operator[]( unsigned p_i ) const;

  /*
   * Function to resize the bucket
   @param size - The new size of the bucket
  */
  void Resize( size_t p_size );   
  
  /*
   * Function to set the pointer to the first object that the container will point to
   @param begin - Pointer to the first object that the container will point to
  */
  void SetStartingIndex( const T* p_begin );

  /*
   * Function that sets the vector of objects to be referenced by this bucket
   * This is used if the objects to be referenced are to be stored only temporarily
   @param values - The vector holding the objects to be referenced
  */
  void SetValues( const std::vector<T>& p_values );
};

    
template<typename T>
Bucket<T>::Bucket( size_t p_capacity ) : _size( 0 ), _capacity( p_capacity ), _beginData( nullptr ) {
}


template<typename T>
size_t
Bucket<T>::GetSize() const {
  return _size;
}


template<typename T>
size_t
Bucket<T>::GetCapacity() const {
  return _capacity;
}


template<typename T>
const T&
Bucket<T>::operator[]( unsigned p_i ) const {
  return *( _beginData + p_i );
}
    

template<typename T>
void
Bucket<T>::Resize( size_t p_size ) {
  _size = p_size;
}
    

template<typename T>
void
Bucket<T>::SetStartingIndex( const T* p_begin ) {
  _beginData = p_begin;
}


template<typename T>
void
Bucket<T>::SetValues( const std::vector<T>& p_values ) {
  _values = p_values;
  _beginData = &_values[0];
}
    
} // Entity
} // Mesh

#endif
