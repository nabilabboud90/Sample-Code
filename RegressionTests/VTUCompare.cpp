#include <unistd.h>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <memory>
#include <vector>
#include <map>

struct Field {
public:
  explicit Field( std::string p_name, int p_components, std::vector<std::string> p_sourceFiles, std::vector<std::string> p_outputType ) : _name( p_name ),
																	  _components( p_components ),
																	  _sourceFiles( p_sourceFiles ),
																	  _outputType( p_outputType ) {}

  void Print() const {
    std::cout << "Field name: " << _name << std::endl;
    std::cout << "Number of components: " << _components << std::endl;
    for( size_t i=0; i<_sourceFiles.size(); ++i ) {
      std::cout << "Source file and output type: " << _sourceFiles[i] << "/" << _outputType[i] << std::endl;
    }
  }

  bool operator==( const Field& p_field ) const {
    if( _name != p_field._name ) { return false; }
    if( _components != p_field._components ) { return false; }

    for( const auto& data : _data ) {
      const auto& coords = data.first;
      const auto& values = data.second;

      if( p_field._data.find( coords ) != p_field._data.end() ) {
	if( p_field._data.find( coords )->second != values ) {
	  return false;
	}
      }
      else {
	return false;
      }
    }

    return true;
  }
  
  std::string _name = {};
  int _components = {};
  std::vector<std::string> _sourceFiles = {};
  std::vector<std::string> _outputType = {};
  std::map<std::vector<double>, std::vector<double> > _data = {};
};


// Function that checks if a file exists
bool IsFileExisting ( const std::string& p_filePath, const std::string& p_name ) {
  if ( FILE *file = fopen( ( p_filePath + "/" + p_name ).c_str(), "r" ) ) {
    fclose(file);
    return true;
  }
  else {
    return false;
  }   
}


std::vector<std::string> ReadPVDFile( const std::string& p_filePath, const std::string& p_fileName ) {
  std::ifstream infile( p_filePath + "/" + p_fileName );

  // Read the pvd file line by line
  std::vector<std::string> pvtu_files_name;
  std::string line;
  while( std::getline( infile, line ) ) {
    // Find the first occurence of the keyword "file="
    size_t found = line.find( "file=" );
    while( found != std::string::npos ) {
      // Check that the next character is a double quotation
      if( line[ found + 5 ] != '"' ) {
	throw std::runtime_error( "The reader expects a \"\"" );
      }

      // Find the closing quotation
      size_t find_quotation = line.find( "\"", found + 6 );

      //Read the name of the pvtu file
      std::string pvtu_file_name( line.begin() + found + 6, line.begin() + find_quotation );
      pvtu_files_name.push_back( pvtu_file_name );
      
      // Find the next occurence of "file=" if any
      found = line.find( "file=", found + 1 );
    }
  }

  return pvtu_files_name;
}


void ReadPVTUFile( const std::string& p_filePath, const std::string& p_fileName, int p_fileIndex, std::map<std::pair<std::string, int>, std::unique_ptr<Field> >& p_data ) {
  std::ifstream infile( p_filePath + "/" + p_fileName );

  // Read the whole file and store it in "file"
  std::string file( "" );
  if( infile ) {
    std::ostringstream ss;
    ss << infile.rdbuf();
    file = ss.str();
  }
  
  // Read the name of the field
  std::vector<std::string> keywords = { "PPointData Scalars=\"", "PCellData Scalars=\"" };
  std::vector<std::string> output_type = { "Point", "Cell" };
  std::string keyword = keywords[0];
  size_t find_field_name_begining = file.find( keyword );
  size_t cnt = 1;
  while( find_field_name_begining == std::string::npos ) {
    keyword = keywords[cnt];
    find_field_name_begining = file.find( keyword );
    cnt++;

    if( cnt > keywords.size() ) {
      throw std::runtime_error( "Trying to read a file from the pvtu file that was written on an entity type that is not supported" );
    }
  }
  
  size_t find_field_name_end = file.find( "\"", find_field_name_begining + keyword.size() );
  std::string field_name( file.begin() + find_field_name_begining + keyword.size(), file.begin() + find_field_name_end );

  // Read the number of components of the field
  keyword = "NumberOfComponents=\"";
  size_t find_num_of_components_begining = file.find( keyword, find_field_name_end + 1 );
  size_t find_num_of_components_end = file.find( "\"", find_num_of_components_begining + keyword.size() );
  std::string num_of_components( file.begin() + find_num_of_components_begining + keyword.size(), file.begin() + find_num_of_components_end );
  
  // Read the name of the file containing the data
  std::vector<std::string> source_files;
  std::vector<std::string> output_types;
  keyword = "Piece Source =\"";
  size_t find_source_file_begining = file.find( keyword, find_num_of_components_end + 1 );
  while( find_source_file_begining != std::string::npos ) {
    size_t find_source_file_end = file.find( "\"", find_source_file_begining + keyword.size() );
    std::string source_file( file.begin() + find_source_file_begining + keyword.size(), file.begin() + find_source_file_end );
    source_files.push_back( source_file );

    // Read the name of another source file if any (this is needed if the problem is ran with multiple mpi processes)
    find_source_file_begining = file.find( keyword, find_source_file_end + 1 );

    // Store the type of the entity on which the field is written for this vtu file
    output_types.push_back( output_type[ cnt - 1 ] );
  }

  // Create and store the field for the corresponding time step
  auto field = std::make_unique<Field>( field_name, std::stoi( num_of_components ), source_files, output_types );
  p_data.insert( std::make_pair( std::make_pair( field_name, p_fileIndex ), std::move( field ) ) );
}


void ReadVTUFile( const std::string& p_filePath, const std::string& p_fileName, Field& p_field, const std::string& p_outputType ) {
  std::ifstream infile( p_filePath + "/" + p_fileName );
  
  // Get some information about the field
  const auto& field_name = p_field._name;
  const auto& num_components = p_field._components;
  
  // Read the size of field data to be read from the current file
  std::string line;
  std::getline( infile, line );
  std::getline( infile, line );
  std::getline( infile, line );
  
  // Read the number of points written in this vtu file
  std::string keyword( "Piece NumberOfPoints = \"" );
  size_t find_num_of_points_begining = line.find( keyword );
  size_t find_num_of_points_end = line.find( "\"", find_num_of_points_begining + keyword.size() );
  int num_of_points = std::stoi( std::string( line.begin() + find_num_of_points_begining + keyword.size(), line.begin() + find_num_of_points_end ) );
  
  // Read and set the number of entities on which the data to be read is defined
  size_t find_num_of_entities_begining, find_num_of_entities_end;
  int num_of_entities = 0;
  if( p_outputType == "Point" ) {
    num_of_entities = num_of_points;
  }
  else if( p_outputType == "Cell" ) {
    keyword = "NumberOfCells = \"";
    find_num_of_entities_begining = line.find( keyword );
    find_num_of_entities_end = line.find( "\"", find_num_of_entities_begining + keyword.size() );
    num_of_entities = std::stoi( std::string( line.begin() + find_num_of_entities_begining + keyword.size(), line.begin() + find_num_of_entities_end ) );
  }
  else {
    throw std::invalid_argument( "The entity type on which the field data is written in the VTU file is not supported by the ReadVTUFile function" );
  }
  
  // Skip the next two lines
  std::getline( infile, line );
  std::getline( infile, line );
  
  // Read the data in the vtu file
  std::vector<std::vector<double> > data;
  int cnt = 0;
  for( unsigned n=0; n<num_of_entities; ++n ) {
    std::vector<double> val;
    
    for( unsigned c=0; c<num_components; ++c ) {
      std::getline( infile, line );
      val.push_back( std::stod( line ) );
    }
    
    data.push_back( val );
  }
  
  // Skip the next four lines
  std::getline( infile, line );
  std::getline( infile, line );
  std::getline( infile, line );
  std::getline( infile, line );

  // Read the coordinates and store them
  std::vector<std::vector<double> > all_points_coordinates;
  for( int c=0; c<num_of_points; ++c ) {
    std::getline( infile, line );

    std::vector<double> coords;
    auto str_stream = std::stringstream( line );
    auto str_iter = std::istream_iterator<double>( str_stream );
    const auto end_of_str = std::istream_iterator<double>();
    auto back_insert_iter = std::back_inserter( coords );
    std::copy( str_iter, end_of_str, back_insert_iter );
    all_points_coordinates.push_back( coords );
  }

  // Store the data in the field
  // If the field is defined at the nodes store the values of the field at each node
  // If the field is stored at the element level, form a vector of vectors where the inner vector holds the coordinates of the points that make up the element and store the values of the field for those points
  if( p_outputType == "Point" ) {
    for( int c=0; c<num_of_points; ++c ) {
      p_field._data.insert( std::make_pair( all_points_coordinates[c], data[ c ] ) );
    }
  }
  else if( p_outputType == "Cell" ) {
    // Skip the next four lines to get to the connectivity section
    std::getline( infile, line );
    std::getline( infile, line );
    std::getline( infile, line );
    std::getline( infile, line );

    for( int c=0; c<num_of_entities; ++c ) {
      std::getline( infile, line );
      
      std::vector<int> connectivity;
      auto str_stream = std::stringstream( line );
      auto str_iter = std::istream_iterator<int>( str_stream );
      const auto end_of_str = std::istream_iterator<int>();
      auto back_insert_iter = std::back_inserter( connectivity );
      std::vector<double> cell_vertex_coords;
      for( int n=0; n<connectivity.size(); ++n ) {
	const auto& vertex_coords = all_points_coordinates[ connectivity[n] ];
	cell_vertex_coords.insert( cell_vertex_coords.end(), vertex_coords.begin(), vertex_coords.end() ); 
      }
      p_field._data.insert( std::make_pair( cell_vertex_coords, data[ c ] ) );
    }
  }
  else {
    throw std::invalid_argument( "The entity type on which the field data is written in the VTU file is not supported by the ReadVTUFile function" );
  }
}


void ReadFile( const std::string& p_filePath, const std::string& p_file, std::map<std::pair<std::string, int>, std::unique_ptr<Field> >& p_dataFile ) {
  // Read the pvd file
  std::vector<std::string> pvtu_files_name = ReadPVDFile( p_filePath, p_file );
  
  // Read the pvtu files
  int cnt = 0;
  for( const auto& name : pvtu_files_name ) {
    ReadPVTUFile( p_filePath, name, cnt, p_dataFile );
    cnt++;
  }

  // Read the vtu files and store the fields data
  for( auto& info : p_dataFile ) {
    const auto& vtu_files = info.second->_sourceFiles;

    unsigned vtu_files_cnt = 0;
    for( const auto& file : vtu_files ) {
      ReadVTUFile( p_filePath, file, *info.second, info.second->_outputType[ vtu_files_cnt ] );
      vtu_files_cnt++;
    }
  }  
}


int main( int argc, char* argv[] ) {
  // Check that the name of 2 files and the path to the directory where they are stored are exactly given to the program
  if( ( argc - 1 ) != 4 ) {
    throw std::invalid_argument( "The comparison function should take exactly the name of two files along with the path to each of the files" );
  }
  
  // Get the names of the files to compare
  std::string file1Path( argv[1] );
  std::string file1( argv[2] );
  std::string file2Path( argv[3] );
  std::string file2( argv[4] );
  
  // Check that both files exist
  if( IsFileExisting( file1Path, file1 ) == false ) {
    throw std::runtime_error( "The file \"" + file1 + " in " + file1Path + "\" does not exist" );
  }
  if( IsFileExisting( file2Path, file2 ) == false ) {
    throw std::runtime_error( "The file \"" + file2 + " in " + file2Path + "\" does not exist" );
  }
  
  // Read and store the data from the first file
  std::map<std::pair<std::string, int>, std::unique_ptr<Field> > data1;
  ReadFile( file1Path, file1, data1 );
  
  // Read and store the data from the second file
  std::map<std::pair<std::string, int>, std::unique_ptr<Field> > data2;
  ReadFile( file2Path, file2, data2 );
  
  // Compare the data in the two files
  std::ofstream output_file( "compare_" + file1 + "_" + file2 + ".log" );
  output_file << "FILE1: " << file1Path + file1 << std::endl;
  output_file << "FILE2 " << file2Path + file2 << std::endl;
  output_file << std::endl;
  output_file << "******************************* START COMPARISON *******************************" << std::endl;
  
  // Check that the two files have the same number of steps
  if( data1.size() != data2.size() ) {
    output_file << "Error! " << file1Path + file1 << " size: " << data1.size() << "   " << file2Path + file2 << " size: " << data2.size() << std::endl;
    output_file.close();
    return 0;
  }
  
  for( const auto& field1_info : data1 ) {
    // Get the pair holding the name of the field and the time step
    const auto& name_timestep = field1_info.first;

    // Get the pointer to the field1
    const auto& field1 = field1_info.second;

    // Get the name of the field being compared
    const auto& field_name = field1_info.first.first;

    // Get the step number at which the two fields will be compared
    const auto& step_num = field1_info.first.second;
    
    // Get the pointer to the field in file2
    if( data2.find( name_timestep ) != data2.end() ) {
      const auto& field2 = data2.find( name_timestep )->second;
      bool are_equal = (*field1 == *field2);
      std::string same = ( are_equal ) ? "true" : "false";
      if( are_equal == true ) {
	output_file << "Field: " << field_name << "   " << "Step: " << step_num << "   " << "Same: " << same << std::endl;
      }
      else {
	output_file << "Error! " << "Field: " << field_name << "   " << "Step: " << step_num << "   " << "Same: " << same << std::endl;
	output_file.close();
	return 0;
      }
    }
    else {
      output_file << "Error! " << "Field: " << field_name << " is not defined in " << file2Path + file2 << std::endl;
      output_file.close();
      return 0;
    }
  }
  output_file << "******************************* END COMPARISON ********************************" << std::endl;
  output_file << "RESULT: THE TWO FILES MATCH" << std::endl;
  
  output_file.close();
  return 1;
}
