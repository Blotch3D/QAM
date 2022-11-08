using Microsoft.Xna.Framework;
using Microsoft.Xna.Framework.Graphics;
using Microsoft.Xna.Framework.Content;
using System;
using Blotch;
using System.Threading;

namespace Blotch.Qam
{
	/// <summary>
	/// The 3D window. This must inherit from BlWindow3D. See BlWindow3D for details.
	/// </summary>
	public class Win3d : BlWindow3D
	{
		public BlSprite TopSprite;
		public BlSprite Sphere;
		public Model SphereModel;
		public SpriteFont Font;

		public static Win3d CreateWindow()
		{
			Win3d win = null;

			new Thread(() =>
			{
				win = new Win3d();
                while (true)
                {
					win.Run();
				}
			}).Start();

			while (win == null)
			{
				Thread.Sleep(10);
			}

			return win;
		}

		/// <summary>
		/// See BlWindow3D for details.
		/// </summary>
		protected override void Setup()
		{
			// Any type of content (3D models, fonts, images, etc.) can be converted to an XNB file by downloading and
			// using the mgcb-editor (see Blotch3D.chm for details). XNB files are then normally added to the project
			// and loaded as shown here. 'Content', here, is the folder that contains the XNB files or subfolders with
			// XNB files. We need to create one ContentManager object for each top-level content folder we'll be loading
			// XNB files from. You can create multiple content managers if content is spread over diverse folders. Some
			// content can also be loaded in its native format using platform specific code (may not be portable) or
			// certain Blotch3D/Monogame methods, like BlGraphicsDeviceManager.LoadFromImageFile.
			var MyContent = new ContentManager(Services, "Content");

			// The model of the toroid
			SphereModel = MyContent.Load<Model>("uv_sphere_192x96");

			Font = MyContent.Load<SpriteFont>("arial14");

			// The sprite we draw in this window
			TopSprite = new BlSprite(Graphics, "Top");
		}

		/// <summary>
		/// See BlWindow3D for details.
		/// </summary>
		/// <param name="timeInfo">Provides a snapshot of timing values.</param>
		protected override void FrameDraw(GameTime timeInfo)
		{
			// Handle the standard mouse and keystroke functions. (This is very configurable)
			Graphics.DoDefaultGui();

			//
			// Draw things here using BlSprite.Draw(), graphics.DrawText(), etc.
			//

			TopSprite.Draw();
		}
	}
}