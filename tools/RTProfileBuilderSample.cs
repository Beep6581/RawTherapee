#region Usings
using System;
using System.Text;
using System.IO;
using System.Globalization;
#endregion

// Raw Therapee sample Custom Profile builder (version 2010-11-06)
//
// How to use:
// 1. In RT, build a profile with the static settings (e.g. Photographers name) and set it as default in the options dialog
// 2. Modify the GetCorrectedSettings function below according to your needs.
// 3. Download and install Microsoft .Net Runtime (latest version is 4.0 as of writing).
//    You can get it for free via Windows Update or from microsoft.com. No need for Visual Studio etc.
// 4. Open a command line and compile this CS-File using the C# 32bit compiler. It is usually installed somewhere here:
//    C:\Windows\Microsoft.NET\Framework\v4.0.30319\csc.exe
//    Call csc.exe with your .CS file as parameter. CSC will compile it and generate an EXE.
// 5. Open your RT options files and find the entry [Profiles]/CustomProfileBuilder.
//    On Win7/Vista it's usually located in C:\Users\<Login>\AppData\Roaming\RawTherapeeAlpha
// 6. Enter the path to your newly built exe here. On Windows, don't forget double slashes (e.g. "C:\\MyDir\\Mybuilder.exe")
//
// And you're done! The EXE is only called on opening the image editor and there is no PP3 yet
// This description is for Windows. The C# code does not use anything fancy, will probably work with MONO on Linux/OSX, too

namespace RTProfilerBuilder {
	/// <summary>Main class. Mostly change GetCorrectedSettings.</summary>
	class RTProfileBuilder {

		/// <summary>This is your personalisation function</summary>
		/// <param name="sectionEntry">Entry, like "Sharpening/Radius"</param>
		/// <param name="value">Current value (from default file)</param>
		/// <param name="fNumber">FNumber</param><param name="exposureSecs">Exposure in seconds</param>
		/// <param name="focalLength">Focal length in MM</param><param name="iso">ISO value</param>
		/// <param name="lens">Lens from EXIF</param><param name="camera">Camera from EXIF</param>
		/// <returns>The value to be written. Simply take the current value if you have nothing to touch.</returns>
		static string GetCorrectedSetting(string sectionEntry, string value,
			float fNumber, float exposureSecs, float focalLength, long iso, string lens, string camera) {

			// We don't do anything to the value if it's not our camera
			if (camera.EndsWith("NIKON D700", StringComparison.InvariantCultureIgnoreCase) && lens.Contains("24.0-120.0")) {
				switch (sectionEntry) {
					// Here is the place to adjust your settings
					// Pretty simple: "SectionName/EntryName" in options file
					// For static parameters it's easier to update the default profile in RT

					case "Vignetting Correction/Amount":
						value = (fNumber < 8 && focalLength < 30) ? "30" : "0";
						break;

					case "RAW/CA":
						value = ToBool(fNumber < 11);  // Means "Enabled if fnumber<11, otherwise disabled"
						break;

					case "Impulse Denoising/Enabled":
						value = ToBool(iso >= 3200);
						break;

					case "HLRecovery/Enabled":
						value = ToBool(iso >= 1600);  // Dynamic range decreases, so we'll probably need it
						break;

					case "Color Boost/Amount":
						if (iso >= 6400) value = "0";  // Colors will get poppy anyway...
						break;

					// Add other parameters here. Mention this is case sensitive!

					default: break;  // we don't touch values we don't care about
				}
			}

			return value;
		}

		#region * Main and Helpers
		static string ToBool(bool condition) { return condition ? "true" : "false"; }
		static string ToFloat(float f) { return f.ToString(CultureInfo.InvariantCulture); }

		/// <summary>Reads default file and parses it. No need to touch it for your personal settings.</summary>
		/// <param name="args">Command line args</param>
		/// <return>0 on all OK.</return>
		static int Main(string[] args) {
			int exitCode = 0;

			try {
				#region Parse input parameters
				int argNo = 0;

				// Name of RAW/JPG to process
				string sourceFile = args[argNo++];

				// What the user selected as his base profile
				string defaultProfileFilePath = args[argNo++];

				// Note that old C++ has no automatic number globalization
				float fNumber = float.Parse(args[argNo++], CultureInfo.InvariantCulture);
				float exposureSecs = float.Parse(args[argNo++], CultureInfo.InvariantCulture);
				float focalLength = float.Parse(args[argNo++], CultureInfo.InvariantCulture);
				long iso = long.Parse(args[argNo++], CultureInfo.InvariantCulture);

				string lens = args[argNo++];
				string camera = args[argNo++];
				#endregion

				// Read default file as basis
				string[] lines = File.ReadAllLines(defaultProfileFilePath);

				// File should be Windows ANSI
				using (TextWriter tw = new StreamWriter(sourceFile + ".pp3", false, Encoding.Default)) {
					string section = "";

					foreach (string line in lines) {
						string l = line.Trim();
						if (!String.IsNullOrWhiteSpace(line)) {

							if (l.StartsWith("["))
								section = l.Trim(new char[] { '[', ']' });
							else if (char.IsLetterOrDigit(l[0]) && l.Contains("=")) {
								int valPos = l.IndexOf("=") + 1;

								string newValue = GetCorrectedSetting(section + "/" + l.Substring(0, valPos - 1), l.Substring(valPos).Trim(),
									fNumber, exposureSecs, focalLength, iso, lens, camera);

								// Merge in new value
								l = l.Substring(0, valPos) + (newValue ?? "");
							}
						}

						tw.WriteLine(l);
					}
				}

			} catch (Exception ex) {
				Console.WriteLine("Error: " + ex.ToString());  // can be seen in the RT console window

				exitCode = 1;
			}

			return exitCode;
		}
		#endregion
	}
}
