# -*- encoding: utf-8 -*-
$:.push File.expand_path("../lib", __FILE__)
require "osgb_convert/version"

Gem::Specification.new do |s|
  s.name        = "osgb_convert"
  s.version     = OsgbConvert::VERSION
  s.platform    = Gem::Platform::RUBY
  s.authors     = ["Matt Patterson"]
  s.email       = ["matt@reprocessed.org"]
  s.homepage    = "http://github.com/fidothe/osgb_convert"
  s.summary     = %q{Geo coordinate transformation between WGS84 (GPS) and OSGB36 (UK Ordnance Survey mapping)}
  s.description = %q{Provides a simple interface to transform Geographic coordinates between WGS84 (GPS) and OSGB36 (UK Ordnance Survey mapping)}

  s.files         = `git ls-files`.split("\n")
  s.test_files    = `git ls-files -- {test,spec,features}/*`.split("\n")
  s.executables   = `git ls-files -- bin/*`.split("\n").map{ |f| File.basename(f) }
  s.require_paths = ["lib"]
end
