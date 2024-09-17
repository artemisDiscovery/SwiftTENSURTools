
import XCTest
@testable import SwiftTENSURTools

import Foundation 
import MathTools

enum tensurError: Error {
    case testImportError
    case parseError
}


func readCoordinates( _ lines:[String], _ index : inout Int ) throws -> [Vector] {
    let ncoord = Int(lines[index])!

    var coords = [Vector]()


        for _ in 0..<ncoord {
            index += 1
            let tokens = lines[index].split(separator: "\t") 
                .map { $0.trimmingCharacters(in: .whitespacesAndNewlines)}
            //print("tokens : \(tokens)")
            
            let x = try Double(tokens[0])
            let y = try Double(tokens[1])
            let z = try Double(tokens[2])

            if x == nil || y == nil || z == nil {
                throw tensurError.parseError
            }

            let coord = Vector([x!,y!,z!])
            coords.append(coord)
            
            
        }
    
    // get to next line 
    index += 1
    return coords 
}

func readRadii( _ lines:[String], _ index : inout Int ) -> [Double] {
    let nrad = Int(lines[index])!

    var radii = [Double]()

    for _ in 0..<nrad {
        index += 1
        let rad = Double(lines[index].trimmingCharacters(in: .whitespacesAndNewlines))!
        radii.append(rad)
    }
    // get to next line 
    index += 1
    return radii 
}

func readCircles(_ lines:[String], _ index : inout Int) -> ([Double],[[AtomCircle]]) {
    let tokens = lines[index].split(separator:"\t")
    let nlevels = Int(tokens[0])!

    var levels = [Double]()
    var circlesForLevels = [[AtomCircle]]()

    for _ in 0..<nlevels {
        index += 1
        let tokens = lines[index].split(separator:"\t")
            .map { $0.trimmingCharacters(in: .whitespacesAndNewlines)}
        let level = Double(tokens[0])!
        levels.append(level)
        var circles = [AtomCircle]()
        let ncircles = Int(tokens[1])!

        for _ in 0..<ncircles {
            index += 1
            var tokens = lines[index].split(separator:"\t")
                .map { $0.trimmingCharacters(in: .whitespacesAndNewlines)}
            let atom = Int(tokens[0])!
            index += 1
            tokens = lines[index].split(separator:"\t")
                .map { $0.trimmingCharacters(in: .whitespacesAndNewlines)}
            let x = Double(tokens[0])!
            let y = Double(tokens[1])!
            let z = Double(tokens[2])!
            let r = Double(tokens[3])!
            
            circles.append(AtomCircle(atom, Vector([x,y,z]), r, AXES.X))
        }

        circlesForLevels.append(circles)
    }

    index += 1 

    return (levels,circlesForLevels)

}

func readContours(_ lines:[String], _ index : inout Int) -> [[[(Int,Vector?,Int?,Vector?,Int?)]]] {
    let tokens = lines[index].split(separator:"\t")
        .map { $0.trimmingCharacters(in: .whitespacesAndNewlines)}
    let nlevels = Int(tokens[0])!

    var contoursForLevels = [[[(Int,Vector?,Int?,Vector?,Int?)]]]()

    for _ in 0..<nlevels {
        index += 1
        var tokens = lines[index].split(separator:"\t")
            .map { $0.trimmingCharacters(in: .whitespacesAndNewlines)}
        
        // contour = [(circleidx,Vector?,Int?,Vector?,Int?)]
        var contours = [[(Int,Vector?,Int?,Vector?,Int?)]]()
        let ncontours = Int(tokens[1])!

        for _ in 0..<ncontours {
            index += 1
            tokens = lines[index].split(separator:"\t")
                .map { $0.trimmingCharacters(in: .whitespacesAndNewlines)}
            let narcs = Int(tokens[0])!
            var contour = [(Int,Vector?,Int?,Vector?,Int?)]()
            for _ in 0..<narcs {
                index += 1
                tokens = lines[index].split(separator:"\t")
                    .map { $0.trimmingCharacters(in: .whitespacesAndNewlines)}
                let cidx = Int(tokens[0])!
                if lines[index].contains("nonsingleton") {
                    index += 1
                    tokens = lines[index].split(separator:"\t")
                        .map { $0.trimmingCharacters(in: .whitespacesAndNewlines)}
                    let xI = Double(tokens[0])!
                    let yI = Double(tokens[1])!
                    let zI = Double(tokens[2])!
                    let circleI = Int(tokens[3])!

                    index += 1
                    tokens = lines[index].split(separator:"\t")
                        .map { $0.trimmingCharacters(in: .whitespacesAndNewlines)}
                    let xJ = Double(tokens[0])!
                    let yJ = Double(tokens[1])!
                    let zJ = Double(tokens[2])!
                    let circleJ = Int(tokens[3])!

                    contour.append( (cidx,Vector([xI,yI,zI]),circleI,Vector([xJ,yJ,zJ]),circleJ) )
                }
                else {
                    contour.append( (cidx,nil,nil,nil,nil) )
                }
            }
            contours.append(contour)
        }

        contoursForLevels.append(contours)
    }

    index += 1 

    return contoursForLevels 

}


func importData(_ path:String) throws -> 
    ([Vector]?,[Double]?,Double?,[Double]?,[[AtomCircle]]?,[[[(Int,Vector?,Int?,Vector?,Int?)]]]?) {

    var input:[String]?

    do {
        input = try String(contentsOfFile:path).split { $0 .isNewline } .map { String($0) }
    }
    catch {
        print("exception reading data file")
        throw tensurError.testImportError
    }
    
    let lines = input!


    var lidx = 0 
    var coordinates:[Vector]?
    var radii:[Double]?
    var probeRad:Double?
    var circlesForLevels:[[AtomCircle]]?
    var levels:[Double]?

    var contoursForLevels:[[[(Int,Vector?,Int?,Vector?,Int?)]]]?

    while lidx < lines.count {
        if lines[lidx].contains("# coordinates") {
            lidx += 1
            do {
                coordinates = try readCoordinates(lines, &lidx )
            }
            catch {
                print("exception readCoordinates")
            }
            
            continue
        }

        if lines[lidx].contains("# radii") {
            lidx += 1
            radii = readRadii(lines, &lidx )
            continue
        }

        if lines[lidx].contains("# probeRad") {
            lidx += 1
            probeRad = Double(lines[lidx].trimmingCharacters(in: .whitespacesAndNewlines))
            lidx += 1
            continue
        }

        if lines[lidx].contains("# atom circles") {
            lidx += 1
            let cdata = readCircles(lines, &lidx )
            levels = cdata.0 
            circlesForLevels = cdata.1
            continue
        }

        if lines[lidx].contains("# contours") {
            lidx += 1
            contoursForLevels = readContours(lines, &lidx )
            
            continue
        }

        lidx += 1


    }

    return (coordinates,radii,probeRad,levels,circlesForLevels,contoursForLevels)


}

final class SwiftTENSURToolsTests: XCTestCase {
    func testExample() throws {
        // XCTest Documentation
        // https://developer.apple.com/documentation/xctest

        // Defining Test Cases and Test Methods
        // https://developer.apple.com/documentation/xctest/defining_test_cases_and_test_methods
    }

    func testCirclesArcs() {

        // make some circles with simple geometry, test intersections

        let PI = acos(-1.0)

        var axis = AXES.X 

        let centerA = Vector([0.0,0.0,0.0])
        let centerB = Vector([0.0,2.0,0.0])
        let centerC = Vector([0.0,1.0,2.0*sin(60.0*(PI/180.0))])

        var circleA = AtomCircle(0, centerA, 1.25,  axis )
        var circleB = AtomCircle(1, centerB, 1.25,  axis )
        var circleC = AtomCircle(2, centerC, 1.25,  axis )

        print("Before intersection : ")
        print("A: \(circleA.str())")
        print("B: \(circleB.str())")
        print("C: \(circleC.str())")

        // intersect A and B

        intersectAtomCircles( circleA, circleB )

        print("\nAfter intersecting A and B : ")

        print("A: \(circleA.str())")
        print("B: \(circleB.str())")

        // intersect C with A and B

        intersectAtomCircles( circleA, circleC )
        intersectAtomCircles( circleB, circleC )

        print("\nAfter second intersection, A and B with C  : ")

        print("A: \(circleA.str())")
        print("B: \(circleB.str())")
        print("C: \(circleC.str())")

        // check whether vectors inside arcs

        let uAC = centerC.sub(centerA).unit()!
        var inside = circleA.exposure[0].contains(uAC)
        print("\nuAC inside A exposure = \(inside)")
        inside = circleC.exposure[0].contains(uAC)
        print("\nuAC inside C exposure = \(inside)")
        inside = circleB.exposure[0].contains(uAC)
        print("\nuAC inside B exposure = \(inside)")

        // add small circle D at upper intersection of A and B, center =  0.0, 1.0, 0.75
        // should not affect exposures

        let centerD = Vector([0.0, 1.0, 0.75])
        var circleD = AtomCircle(3, centerD, 0.25,  axis )

        intersectAtomCircles( circleA, circleD )
        intersectAtomCircles( circleB, circleD )
        intersectAtomCircles( circleC, circleD )
        
        print("\nAfter adding circle D, intersecting with existing A,B,C  : ")

        print("A: \(circleA.str())")
        print("B: \(circleB.str())")
        print("C: \(circleC.str())")
        print("D: \(circleD.str())")

        
        let centerE = Vector([0.0, 1.0, -0.75])
        var circleE = AtomCircle(4, centerE, 0.25,  axis )

        intersectAtomCircles( circleA, circleE )
        print("After intersecting A and E :")
        print("A: \(circleA.str())")
        print("E: \(circleE.str())")

        intersectAtomCircles( circleB, circleE )
        print("After intersecting B and E :")
        print("B: \(circleB.str())")
        print("E: \(circleE.str())")

        //
        intersectAtomCircles( circleC, circleE )
        print("After intersecting C and E :")
        print("C: \(circleC.str())")
        print("E: \(circleE.str())")
        intersectAtomCircles( circleD, circleE )
        
        print("\nAfter adding circle E, intersecting with existing A,B,C,D  : ")

        print("A: \(circleA.str())")
        print("B: \(circleB.str())")
        print("C: \(circleC.str())")
        print("D: \(circleD.str())")    
        print("E: \(circleE.str())")

        var circles = [circleA,circleB,circleC,circleE]

        var contour:Contour?

        do {
            contour = try Contour(circles)
        }
        catch {
            print("exception in Contour")
        }

        print(contour!.str(verbose:true))

        // Regenerate circles, try pairwise

        print("pairwise intersection")

        circleA = AtomCircle(0, centerA, 1.25,  axis )
        circleB = AtomCircle(1, centerB, 1.25,  axis )
        circleC = AtomCircle(2, centerC, 1.25,  axis )
        circleD = AtomCircle(3, centerD, 0.25,  axis )
        circleE = AtomCircle(4, centerE, 0.25,  axis )

        circles = [circleA,circleB,circleC,circleD,circleE]

        for i in 0..<(circles.count-1) {
            for j in (i+1)..<circles.count {
                intersectAtomCircles(circles[i],circles[j])
            }
        }

        do {
            contour = try Contour(circles)
        }
        catch {
            print("exception in Contour")
        }

        print(contour!.str(verbose:true))

        print("pairwise intersection, change order")

        circleA = AtomCircle(0, centerA, 1.25,  axis )
        circleB = AtomCircle(1, centerB, 1.25,  axis )
        circleC = AtomCircle(2, centerC, 1.25,  axis )
        circleD = AtomCircle(3, centerD, 0.25,  axis )
        circleE = AtomCircle(4, centerE, 0.25,  axis )

        circles = [circleC,circleE,circleA,circleD,circleB]

        for i in 0..<(circles.count-1) {
            for j in (i+1)..<circles.count {
                intersectAtomCircles(circles[i],circles[j])
            }
        }

        do {
            contour = try Contour(circles)
        }
        catch {
            print("exception in Contour")
        }

        print(contour!.str(verbose:true))

    }

    
    
    
    func testLayers() throws {

        


        // Make collection of atom circles

        
        // read circle definitions from data directory
        
        var importdata:([Vector]?,[Double]?,Double?,[Double]?,[[AtomCircle]]?,[[[(Int,Vector?,Int?,Vector?,Int?)]]]?)?
        
        do {
            importdata = try importData("Tests/DATA/data.txt")
        }
        catch {
            print("exception calling importData")
            throw tensurError.testImportError
        }

        let data = importdata!

        let coordinates = data.0!
        let radii = data.1!
        let probeRad = data.2!
        let levels = data.3!
        let circlesForLevels = data.4!
        let contoursForLevels = data.5!

        var atomcoords = [Double]() 

        for vec in coordinates {
            atomcoords += vec.coords 
        }

        let atompos = Matrix<Double>([coordinates.count,3], content:atomcoords )

        let maxrad = radii.max()!

        let xcoords = coordinates.map { $0.coords[0] }
        let minX = xcoords.min()!

        let mincoord = levels.min()!

        print("min layer coord this data = \(mincoord)")
        print("min X-coord = \(minX), max rad = \(maxrad), min X-coord - maxrad - probeRad = \(minX - maxrad - probeRad)")

        let circleLayers = atomCirclesForLayers( atompos:atompos, radii:radii, 
            proberad:probeRad, minaxiscoord:mincoord, layerdelta:0.5, axis:AXES.X, numthreads:1 )
            
        // Get all atom circles with X 
        
        var atomcircles = circleLayers.objects as! [AtomCircle]

        let X = 26.5

        var layercircles = atomcircles.filter { $0.center.coords[0] == X }
        
        print("have \(layercircles.count) circles for layer at \(X)")


        //print("layer circle data ---")

        //for circle in layercircles {
        //    print(circle.str())
        //}

        for i in 0..<(layercircles.count-1) {
            for j in (i+1)..<layercircles.count {
                //let d = layercircles[i].center.sub(layercircles[j].center).length()
                //if d > layercircles[i].radius + layercircles[j].radius {
                //    continue
                //}
                
                    //print("\n\nintersect \(i) and \(j)\n")
                    //print("prior intersect, arcs \(i)---")
                    //for arc in layercircles[i].exposure {
                    //    print(arc.str())
                    //}
                    //print("prior intersect, arcs \(j)---")
                    //for arc in layercircles[j].exposure {
                    //    print(arc.str())
                    //}
                
                
                intersectAtomCircles(layercircles[i],layercircles[j])
                
                    //print("after intersect, arcs \(i)---")
                    //for arc in layercircles[i].exposure {
                    //    print(arc.str())
                    //}
                    //print("after intersect, arcs \(j)---")
                    //for arc in layercircles[j].exposure {
                    //    print(arc.str())
                    //}
                
            }
        }

        // look for singleton circle

        let UP = axisUP[AXES.X.rawValue]
        let RIGHT = axisRIGHT[AXES.X.rawValue]

        //print("singleton circle")

        for circle in layercircles {
            if circle.center.coords[UP] > 5.0 && circle.center.coords[RIGHT] < 18.0 {
                print(circle.str())
                //print("exposure: \(circle.exposure.count)")
            }
        }

        var contour:Contour?

        do {
            contour = try Contour(layercircles)
        }
        catch {
            print("exception in Contour first round")
        }

        if contour != nil {
            //print("CONTOUR, first round ----- ")
            //print(contour!.str(verbose:true))
        }
        else {
            print("contour failed first round")
        }

        //print("\ncontour python ...\n")

        //let text = printPython( layercircles, contour!.arcsInOrder)

        //print(text)

        // check for removed arcs - OK working now

        
        //print("removed circles/arcs ---")

        for circle in layercircles {
            if circle.removed {
                //print("removed circle at atom \(circle.atom)")
                continue
            }
            else {
                //print("retained circle at atom \(circle.atom)")
            }
            for arc in circle.exposure {
                if arc.removed {
                    //print("removed arc")
                    //print(arc.str())
                }
            }
        }
       

        // print all arcs 

        //print("all arcs :")

        //for circle in layercircles {
        //    if circle.removed {
        //        continue 
        //    }
        //    for arc in circle.exposure {
        //        print(arc.str())
        //    }
        //}
        
        

        //var axis:AXES

        //var mincoord:Double
        //var maxcoord:Double

        //var delta:Double

        //var objects:[Any]

        //var layers:[Int]

        //var layerBits:[objectBits?]

        print("layers min, max coord = \(circleLayers.mincoord), \(circleLayers.maxcoord)")
        print("layers count = \(circleLayers.layerBits.count)")
        print("have \(circleLayers.objects.count) atom circles")
        print("\n\n")

        
 
        // Get all atomcircles inside another and remove 
    

        atomcircles = circleLayers.objects as! [AtomCircle]
        var centers = [Double]() 
        var centerradii = atomcircles.map { $0.radius }

        for atomcircle in atomcircles {
            centers += atomcircle.center.coords
        }

        print("have \(atomcircles.count) atom circles, length centers array = \(centers.count)")

        for lidx in 0..<circleLayers.layerBits.count {
            let lb = circleLayers.layerBits[lidx]
            if lb == nil {
                //print("layer \(lidx) : nil ")
            }
            else {
                let indices = lb!.indices()
                let first = indices[0]
                let x = atomcircles[first].center.coords[0]
                let y = atomcircles[first].center.coords[1]
                let z = atomcircles[first].center.coords[2]
                //print("layer \(lidx) : \(lb!.indices().count) atom circles, center at \(x), \(y), \(z) ")
            }
            
        }

        
        let surfdata = generateSurfaceProbes( coordinates:coordinates, radii:radii, probeRadius:probeRad, levelspacing:0.5, minoverlap:0.5, numthreads:10,
                         skipCCWContours:true )

        let probes = surfdata.0 

        var str = ""
        

        for probe in probes {
            if !probe.singleton {
                str += "\(probe.center.coords[0]) \(probe.center.coords[1]) \(probe.center.coords[2]) \(probeRad) 1.0 1.0 1.0\n"
            }
            else {
                str += "\(probe.center.coords[0]) \(probe.center.coords[1]) \(probe.center.coords[2]) \(probeRad) 1.0 0.0 0.0\n"
            }
        }

        do {
            try str.write(to: URL(fileURLWithPath:"probes.txt"), atomically: true, encoding: String.Encoding.utf8)
        } catch {
            print("could not write probes file !")
        }
        
        let time1 = Date().timeIntervalSince1970

        var tridata:([Vector],[Vector],[[Int]])?

        // set axis for testing 

        let testaxis = AXES.Z

        do {
            tridata = try generateTriangulation( probes:probes, probeRadius:probeRad, gridspacing:0.4, 
            densityDelta:0.1, densityEpsilon:0.1, isoLevel:1.0, numthreads:10, axis:testaxis) 
        }
        catch {
            print("triangulation code failed !")
        }

    // write out in obj format

    let VERTICES = tridata!.0 
    let NORMALS = tridata!.1
    let FACES = tridata!.2

    let url = URL(fileURLWithPath: "./isosurface.axis.\(testaxis.rawValue).obj")

    var outstr = ""

    for vertex in VERTICES {
        outstr += "v \(vertex.coords[0]) \(vertex.coords[1]) \(vertex.coords[2])\n"
    }

    for normal in NORMALS {
        outstr += "v \(normal.coords[0]) \(normal.coords[1]) \(normal.coords[2])\n"
    }
    
    for face in FACES {
        outstr += "f \(face[0]+1) \(face[1]+1) \(face[2]+1)\n"
    }
    

    do {
        try outstr.write(to: url, atomically: true, encoding: String.Encoding.utf8)
    } catch {
        print("error writing file \(url)")
    }

    print("wrote file to path \(url)")

    let time2 = Date().timeIntervalSince1970

    print("\ntime for marching cubes = \(time2 - time1)")


        
    }
    
   
    

    

 
}
