 #define itkGetConstStdVectorMacro(name,type) \
  virtual const type & Get##name () \
  { \
    itkDebugMacro("returning " << #name ); \
    return this->m_##name; \
  } 

#define itkSetConstStdVectorMacro(name,type) \
  virtual void Set##name (const type & arg) \
  { \
    itkDebugMacro("setting " << #name ); \
    this->m_##name = arg; \
    this->Modified(); \
  }

#define itkAddToStdVectorMacro(name,type,vector) \
  virtual void Add##name ( const type & arg ) \
  { \
    itkDebugMacro("adding " << arg << " to " << #vector ); \
    this->m_##vector.push_back( arg ); \
    this->Modified(); \
  }

#define itkAddConstObjectToStdVectorMacro(name,type,vector) \
  virtual void Add##name ( const type* arg ) \
  { \
    itkDebugMacro("adding " << arg << " to " << #vector ); \
    this->m_##vector.push_back( arg ); \
    this->Modified(); \
  } 

#define itkAddObjectToStdVectorMacro(name,type,vector) \
  virtual void Add##name ( type* arg ) \
  { \
    itkDebugMacro("adding " << arg << " to " << #vector ); \
    this->m_##vector.push_back( arg ); \
    this->Modified(); \
  } 
