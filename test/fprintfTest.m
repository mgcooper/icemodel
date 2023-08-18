classdef fprintfTest < matlab.perftest.TestCase
   % classdef fprintfTest < matlab.unittest.TestCase

   % This example shows how to create and run a class-based performance test and
   % regression test for the fprintf function.

   % This is Mathworks example for class-based perftest

   % If the test inherits from matlab.unittest.TestCase, then the measured time
   % does not include the time to open and close the file or the assertion
   % because these activities take place inside a TestMethodSetup block, and not
   % inside a Test block. However, the measured time includes the time to
   % perform the verifications. Best practice is to measure a more accurate
   % performance boundary. To do that, inherit from matlab.perftest.TestCase and
   % call the startMeasuring and stopMeasuring methods to create a boundary
   % around the fprintf function call. 
   % 
   % When the perftest framework is used, the measured time includes only the
   % call to fprintf, and the testing framework still evaluates the
   % qualifications.

   % See the commented out sections, which would work with the unittest
   % framework, versus uncommented sections, which work with perftest.
   
   

   properties
      file
      fid
   end
   methods(TestMethodSetup)
      function openFile(testCase)
         testCase.file = tempname;
         testCase.fid = fopen(testCase.file,'w');
         testCase.assertNotEqual(testCase.fid,-1,'IO Problem')

         testCase.addTeardown(@delete,testCase.file);
         testCase.addTeardown(@fclose,testCase.fid);
      end
   end

   methods(Test)
      function testPrintingToFile(testCase)
         textToWrite = repmat('abcdef',1,5000000);

         % Use this for the unittest framework
         % fprintf(testCase.fid,'%s',textToWrite);

         % Use this for the perftest framework
         testCase.startMeasuring();
         fprintf(testCase.fid,'%s',textToWrite);
         testCase.stopMeasuring();

         % This is used for both cases
         testCase.verifyEqual(fileread(testCase.file),textToWrite)
      end

      function testBytesToFile(testCase)
         textToWrite = repmat('tests_',1,5000000);

         % Use this for the unittest framework
         % nbytes = fprintf(testCase.fid,'%s',textToWrite);

         % Use this for the perftest framework
         testCase.startMeasuring();
         nbytes = fprintf(testCase.fid,'%s',textToWrite);
         testCase.stopMeasuring();

         % This is used for both cases
         testCase.verifyEqual(nbytes,length(textToWrite))
      end
   end
end