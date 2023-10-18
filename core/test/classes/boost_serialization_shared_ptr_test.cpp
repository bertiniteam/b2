// i developed this code to ask for help on my problem in compilation of SLP's after de-serialization.  namely, bad weak pointers when calling shared_from_this.
//
// i posted a version of this as a MWE to https://stackoverflow.com/questions/77068663/bad-weak-pointer-after-deserializing-object-that-has-a-shared-pointer-to-base
// user @sehe provided an answer: move the inheritance from `enable_shared_from_this` to `Base` from the derived type. 
//
// thus, i think my problem with SLP's can be solved by moving the inheritance

#include <memory>


#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/nvp.hpp>

//this #define MUST appear before #include <boost/test/unit_test.hpp>
// #define BOOST_TEST_MODULE "MWE" // for standalone compilation
// #define BOOST_TEST_DYN_LINK // for standalone compilation

#include <boost/test/unit_test.hpp>
#include <fstream>








// abstract, and owns a pointer to Base.
class Base: public std::enable_shared_from_this<Base>
{
public:
	Base() = default;
	virtual ~Base() = default;


	void set_ptr(std::shared_ptr<Base> const& n){
		this->_next = n;
	}

	const std::shared_ptr<Base> get_ptr() const{
		return this->_next;
	}

	virtual void do_thing_using_shared_from_this(){
		throw std::runtime_error("failed to override, this is the base");
	}

	virtual void print(){
		std::cout << "this is base, this should have been overridden" << std::endl;
	}

private:

	std::shared_ptr<Base> _next = nullptr;


	friend class boost::serialization::access;

	template <typename Archive>
	void serialize(Archive& ar, const unsigned version) {
		ar & _next;
	}



};

// concrete, with a method using `shared_from_this`
// this is a proxy for a visitable type.
class Node: public /*virtual*/ Base
{
public:
	Node() = default;
	virtual ~Node() = default;
		
	void check_shared_from_this() const{
		std::shared_ptr<const Node> another_shared_ptr = std::static_pointer_cast<Node const>(this->shared_from_this());
	}

	virtual void do_thing_using_shared_from_this() override
	{
		auto as_shared = this->shared_from_this();
	}

	virtual void print() override{
		std::cout << "this is derived" << std::endl;
	}
private:




	friend class boost::serialization::access;

	template <typename Archive>
	void serialize(Archive& ar, const unsigned version) {
		ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Base);
	}

};


BOOST_CLASS_EXPORT(Base)
BOOST_CLASS_EXPORT(Node)


// this is a proxy for a type that owns visitable pointers.
class System{

public:
	System(std::shared_ptr<Base> const& p): _ptr(p) 
	{}

	System() = default;

	~System() = default;

	void set_ptr(std::shared_ptr<Base> const& n){
		this->_ptr = n;
	}

	void check_shared_from_this() const{
		_ptr->do_thing_using_shared_from_this(); // this call fails after de-serialization, due to bad weak pointers
	}

	void call_a_virtual_function(){
		_ptr->print();
	}
private:
	std::shared_ptr<Base> _ptr = nullptr;


	friend class boost::serialization::access;

	template <typename Archive>
	void serialize(Archive& ar, const unsigned version) {
		ar.template register_type<Node>(); // have to register, because have pointer to base
		ar & _ptr; // serialize the data this type owns
	}


};


BOOST_AUTO_TEST_SUITE(boost_serialization)



BOOST_AUTO_TEST_CASE(serialize_deserialize_has_a)
{
	{ // to create a scope


		auto a_ptr = std::make_shared<Node>();
		auto b_ptr = std::make_shared<Node>();

		b_ptr->set_ptr(a_ptr);

		System sys(b_ptr);

		std::ofstream fout("serialization_basic");
		
		boost::archive::text_oarchive oa(fout);
		
		sys.check_shared_from_this();
		// write class instance to archive
		oa << sys;
	}
	
	
	{
		std::ifstream fin("serialization_basic");
		
		boost::archive::text_iarchive ia(fin);
		// read class state from archive

		System sys;
 
		ia >> sys;

		sys.call_a_virtual_function(); // this call is fine
		sys.check_shared_from_this(); // bad weak pointer, due to shared_from_this
	}

}

BOOST_AUTO_TEST_SUITE_END()



// this test doesn't fail, but i think 
BOOST_AUTO_TEST_CASE(serialize_deserialize_shared_ptr_holder)
{
	{ // to create a scope


		auto a_ptr = std::make_shared<Node>();
		auto b_ptr = std::make_shared<Node>();

		b_ptr->set_ptr(a_ptr);

		auto b_check = b_ptr->shared_from_this();
		auto a_check = a_ptr->shared_from_this();

		std::ofstream fout("serialization_basic");
		
		boost::archive::text_oarchive oa(fout);
		
		// write class instance to archive
		oa << b_ptr;
	}
	
	
	{
		std::ifstream fin("serialization_basic");
		
		boost::archive::text_iarchive ia(fin);
		// read class state from archive

		std::shared_ptr<Node> result;

		ia >> result;

		result->check_shared_from_this();
	}

}


