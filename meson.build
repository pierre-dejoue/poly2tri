project('poly2tri', ['cpp'])

include = include_directories('poly2tri/include')
lib = static_library('poly2tri', sources : [
	'poly2tri/src/common/node.cpp',
	'poly2tri/src/common/shapes.cpp',
	'poly2tri/src/sweep/advancing_front.cpp',
	'poly2tri/src/sweep/cdt.cpp',
	'poly2tri/src/sweep/policy.cpp',
	'poly2tri/src/sweep/sweep.cpp',
	'poly2tri/src/sweep/sweep_context.cpp',
], include_directories: include)

thread_dep = dependency('threads')
boost_test_dep = dependency('boost', modules : [ 'filesystem', 'unit_test_framework' ], required : false)
if boost_test_dep.found()
	test('Unit Test', executable('unittest', [
		'unittest/test_cdt.cpp',
		'unittest/test_triangle.cpp',
		'unittest/utility.cpp',
	], dependencies : [boost_test_dep, thread_dep], link_with : lib))
endif

poly2tri_dep = declare_dependency(include_directories : include, link_with : lib)
