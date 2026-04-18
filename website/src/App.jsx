import { useEffect, useState } from 'react';

const releaseBaseUrl = 'https://github.com/PI-Prasaad-Krishna/FormulaOne-replayer/releases/latest';
const issuesUrl = 'https://github.com/PI-Prasaad-Krishna/FormulaOne-replayer/issues';
const feedbackUrl='https://github.com/PI-Prasaad-Krishna/FormulaOne-replayer/discussions';

const downloads = [
  {
    name: 'Windows',
    filename: 'F1 Visualizer Setup.msi',
    details: 'Recommended for Windows 10 and 11. Fast install, native desktop shortcuts, ready to run.',
    badge: 'Most downloaded',
    accent: 'linear-gradient(135deg, rgba(255, 59, 48, 0.24), rgba(255, 255, 255, 0.08))',
    href: releaseBaseUrl,
  },
  {
    name: 'macOS',
    filename: 'F1 Visualizer.dmg',
    details: 'Built for Mac users who want a clean install and smooth launch on Apple Silicon and Intel.',
    badge: 'macOS ready',
    accent: 'linear-gradient(135deg, rgba(255, 255, 255, 0.18), rgba(255, 255, 255, 0.05))',
    href: releaseBaseUrl,
  },
  {
    name: 'Linux',
    filename: 'F1 Visualizer.AppImage',
    details: 'A portable Linux build for quick setup across supported distributions.',
    badge: 'Portable',
    accent: 'linear-gradient(135deg, rgba(0, 194, 255, 0.22), rgba(255, 255, 255, 0.06))',
    href: releaseBaseUrl,
  },
];

const highlights = [
  {
    title: 'Latest release',
    value: 'Available now',
  },
  {
    title: 'Platforms supported',
    value: 'Windows, macOS, Linux',
  },
  {
    title: 'Install experience',
    value: 'Fast and polished',
  },
];

const features = [
  {
    title: 'Immersive Real-Time Telemetry',
    text: 'Experience race data live. Track speed, gears, and braking with stunning visualization widgets.',
  },
  {
    title: 'Cross-Platform Performance',
    text: 'Optimized flawlessly for Windows, macOS, and Linux to deliver a butter-smooth 60fps experience.',
  },
  {
    title: 'Premium Dark-Mode Aesthetics',
    text: 'Built with a sleek, aerodynamic interface featuring dynamic glassmorphism and neon racing accents.',
  },
];

function App() {
  const [mounted, setMounted] = useState(false);

  useEffect(() => {
    setMounted(true);

    const observer = new IntersectionObserver(
      (entries) => {
        entries.forEach((entry) => {
          if (entry.isIntersecting) {
            entry.target.classList.add('is-visible');
          }
        });
      },
      {
        threshold: 0.2,
        rootMargin: '0px 0px -10% 0px',
      },
    );

    document.querySelectorAll('[data-reveal]').forEach((node) => observer.observe(node));

    return () => observer.disconnect();
  }, []);

  return (
    <div className={`site ${mounted ? 'is-ready' : ''}`}>
      <div className="ambient ambient-left" />
      <div className="ambient ambient-right" />
      <div className="grid-overlay" />

      <header className="topbar">
        <a className="brand" href="#top" aria-label="F1 Visualizer home">
          <span className="brand-mark" aria-hidden="true">
            <svg viewBox="0 0 48 48" role="img" aria-hidden="true">
              <path d="M7 27h18l-3 7h9l10-17h-7l4-8H21l-8 10H7z" />
            </svg>
          </span>
          <span>
            <strong>F1 Visualizer</strong>
            <small>Official downloads</small>
          </span>
        </a>

        <nav className="nav">
          <a href="#downloads">Download</a>
          <a href="#features">Features</a>
          <a href={issuesUrl} target="_blank" rel="noreferrer" title="Report issues">Report Issues</a>
          <a href={feedbackUrl} target="_blank" rel="noreferrer" title="Provide feedback">Feedback</a>
        </nav>

        <a className="release-pill" href={releaseBaseUrl} target="_blank" rel="noreferrer">
          Latest release
        </a>
      </header>

      <main id="top">
        <section className="hero">
          <div className="hero-copy" data-reveal>
            <div className="eyebrow">Official release download</div>
            <h1>Download F1 Visualizer for your desktop.</h1>
            <p>
              Get the latest version for Windows, macOS, or Linux. Clean installation, premium
              presentation, and a release built to be easy to trust.
            </p>

            <div className="cta-row">
              <a className="primary-cta" href="#downloads">
                Download latest release
              </a>
              <a className="secondary-cta" href={releaseBaseUrl} target="_blank" rel="noreferrer">
                View on GitHub Releases
              </a>
            </div>

            <div className="stat-row">
              {highlights.map((item) => (
                <div key={item.title} className="stat-card">
                  <span>{item.title}</span>
                  <strong>{item.value}</strong>
                </div>
              ))}
            </div>
          </div>

          <div className="hero-visual" data-reveal>
            <div className="track-panel vehicle-panel">
              <div className="track-label">Formula 1-inspired styling</div>
              <div className="hero-accent-line" aria-hidden="true" />
              <img src="/hero-bg.png" className="hero-image" alt="Formula One style race car illustration" />

              <div className="vehicle-livery" aria-hidden="true">
                <span />
                <span />
                <span />
              </div>

              <div className="build-card">
                <div>
                  <span className="muted">Download destination</span>
                  <strong>GitHub Releases</strong>
                </div>
                <div className="build-meta">
                  <span>Latest version</span>
                  <span>Trusted source</span>
                </div>
              </div>
            </div>
          </div>
        </section>

        <section className="downloads-section" id="downloads">
          <div className="section-heading" data-reveal>
            <span>Downloads</span>
            <h2>Select the installer for your device.</h2>
          </div>

          <div className="download-grid">
            {downloads.map((item) => (
              <a
                key={item.name}
                className="download-card"
                href={item.href}
                target="_blank"
                rel="noreferrer"
                style={{ background: item.accent }}
                data-reveal
              >
                <div className="download-card-top">
                  <span className="download-badge">{item.badge}</span>
                  <span className="download-platform">{item.name}</span>
                </div>

                <div className="download-main">
                  <strong>{item.filename}</strong>
                  <p>{item.details}</p>
                </div>

                <div className="download-footer">
                  <span>Open latest release</span>
                  <span aria-hidden="true">
                    <svg viewBox="0 0 24 24">
                      <path d="M13 5l7 7-7 7M20 12H4" />
                    </svg>
                  </span>
                </div>
              </a>
            ))}
          </div>
        </section>

        <section className="features-panel" id="features" data-reveal>
          <div>
            <span className="eyebrow">Engineered for Performance</span>
            <h2>Everything you need to immerse yourself in the race.</h2>
          </div>

          <div className="step-grid">
            {features.map((feature) => (
              <div key={feature.title} className="step-card">
                <strong>{feature.title}</strong>
                <p>{feature.text}</p>
              </div>
            ))}
          </div>
        </section>
      </main>
    </div>
  );
}

export default App;